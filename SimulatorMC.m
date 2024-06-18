%% Downlink, Monte Carlo of snapshots
% Downlink simulator for Multi-Connectivity.
% The BSs assign a set of RBs to each UE depending on its traffic, and the UEs use all the assigned RBs to transmitt the traffic
% generated following FTP model 3. Multi-Connectivity is performed selecting the UEs with the SNR below a threshold.

clear;
close all;
clc;

%% Parameters
numIterations = 1; % Number of Monte Carlo iterations
Lscenario = 2e3; % Scenario length [m]
scheduling = 'BestCQI'; % Resource scheduling technique, 'RR' or 'BestCQI'
Users = 100; % Number of Users
MCueSNRthresholds = 20; % SNR threshold for MC
NumConn = 2; % Number of MC links
Ptx = 30; % BS transmitted power [dBm]
BLERtarget = 1e-1; % Target BLER
ArrivalRateAvg = 2500; % Average packet arrival rate [pkt/s]
PktSizeAvg = 1e-3; % Average packet size[Mb/pkt]
% BLER curves for the selected MCS
load('BLERcurves3GPP38214.mat');

%% Init
fc = 2; % [GHz]
bandwidth = 20e6; % [Hz]
mu = 0; % 5G numerology.
TTI = 1e-3 / 2^mu; % Transmission Time Interval [s]
numSlotperTTI = 14; % Number of Resource Element per TTI
numSubcarrierPerRB = 12;  % Number of subcarriers per Resource Blocks
K = 100; % Number of available RB per BS
heightUEVar = 21; % UE height variability [m]
heightUEOffset = 1.5; % UE minimum height [m]
heightBS = 25; % BS height variability [m]
heightEff = 1; % Effective height [m]
NoiseSpectralDensity = -174; % Noise spectral density [dBm/Hz]
noisePower = NoiseSpectralDensity + 10*log10(bandwidth); % Noise Power [dBm]
ISD = 500; % [m] Intersite distance. Source 3GPP T.S. 38.901
% Generation of BS with hexagonal cells in the scenario
[bsPosition] = HexGrid([0 0],ISD,Lscenario,heightBS);
bsPosition = [bsPosition(:,1)+Lscenario/2 bsPosition(:,2)+Lscenario/2 bsPosition(:,3)];
bsNumber = size(bsPosition,1);
% Projections of the scenario to consider the boundary effects
Proj = [[0 Lscenario 0];[0 -Lscenario 0];[Lscenario 0 0];[-Lscenario 0 0]; ...
        [Lscenario Lscenario 0];[Lscenario -Lscenario 0]; [-Lscenario Lscenario 0];[-Lscenario -Lscenario 0]];
% Channel model UMa NLOS (Source: 3GPP T.S. 38.901, Table 7.4.1-1)
SFstd = 6; % [dB]
ShadFad = makedist('Normal','mu',0,'sigma',SFstd);
% MCS spectral efficiencies (Source: 3GPP T.S. 38.214, Table 5.1.3.1-2) 
MCSperRE = [0                                                               ... No connection
            0.2344 0.3770 0.6016 0.8770 1.1758                              ... Mod. order 2
            1.4766 1.6953 1.9141 2.1602 2.4063 2.5703                       ... Mod. order 4
            2.7305 3.0293 3.3223 3.6094 3.9023 4.2129 4.5234 4.8164 5.1152  ... Mod. order 6
            5.3320 5.5547 5.8906 6.2266 6.5703 6.9141 7.1602 7.4063         ... Mod. order 8
            ];
% [b/RE] NB: MCSperRE array contains the #bits per RE based on the selected MCS
MCSeff = (MCSperRE.*(numSubcarrierPerRB*numSlotperTTI)*1e-6)/TTI; % Spectral efficiency of each MCS per RB [Mbps/RB]
% MCS thresholds based on SNR and target BLER
SNRthresholds = zeros(1,length(MCSeff));
SNRthresholds(1) = -inf;
for m = 2:length(MCSeff)
    aux = BLERcurves{m-1};
    [idx0,~]=find(aux(:,2)<=BLERtarget,1);
    SNRthresholds(m) = round(aux(idx0,1),3);
end
L = length(MCueSNRthresholds);
% Output variables
result_throughput = cell(numIterations,L);
result_bler = cell(numIterations,L);
result_usage = cell(numIterations,L);
result_MCueList = cell(numIterations,L);

%% Core
disp(strcat('Simulating (',string(datetime),') ...'));
parfor it = 1:numIterations
    %% Init cycle
    bsPos = bsPosition;
    NConn = NumConn;
    Nue = Users;
    Nbs = bsNumber;
    BLcurves = BLERcurves;
    MCS = MCSeff;
    SNRth = SNRthresholds;
    projections = Proj;
    thresholds = MCueSNRthresholds;
    for l = 1:L
        %% Init simulation
        SNRmcTh = thresholds(l);
        availableRatePerUser = zeros(Nue,NConn);
        throughput = zeros(Nue,NConn);
        throughputMC = zeros(Nue,1);
        throughput_aux = 0;
        BLERmc = zeros(Nue,1);
        actualBLERperUser = zeros(Nue,NConn);
        d = single(zeros(Nue,Nbs));
        d2D = single(zeros(Nue,Nbs));
        PL = single(zeros(Nue,Nbs));
        PLlos = single(zeros(Nue,Nbs));
        listRB = zeros(Nbs,2);
        %% UE allocation
        uePos = [rand(Nue,2)*Lscenario, heightUEOffset+rand(Nue,1)*heightUEVar];
        %% Calculate UE-BS distances, PL and SNR
        for u = 1:Nue
            aux1 = sqrt(sum((uePos(u,:)-bsPos(:,:)).^2, 2))';
            aux3 = zeros(1,8);
            for b = 1:Nbs
                % Distance calculation without (if) and with (else) boundary effects
                if aux1(b) <= Lscenario*sqrt(2)/2
                    d(u,b) = single(aux1(b));
                else
                    for i = 1:8
                        aux3(i) = sqrt(sum(((uePos(u,:) + projections(i,:)) - bsPos(b,:)).^2, 2));
                    end
                    d(u,b) = single(min([aux1(b) aux3]));
                end
                % Path loss computation following 3GPP specifications
                dBP = (4*fc*1e9*(heightBS-heightEff)*(uePos(u,3)-heightEff))/physconst('lightspeed');
                if d(u,b)<=dBP % d3D ~= d2D
                    PLlos(u,b) = 28+22*log10(d(u,b))+20*log10(fc);
                else
                    PLlos(u,b) = 28+40*log10(d(u,b))+20*log10(fc)-9*log10((dBP)^2+(heightBS-uePos(u,3))^2);
                end
                PL(u,b) = max([PLlos(u,b) 13.54+39.08*log10(d(u,b))+20*log10(fc)-0.6*(uePos(u,3)-1.5)])+random(ShadFad,1);
            end
        end
        SNR = single(Ptx - PL - noisePower);
    %% MCS assignment for each UE based on its SNR
        CQI = single(zeros(size(SNR)));
        MbpsperRB = single(zeros(size(SNR)));
        for b = 1:Nbs
            for u = 1:Nue
                [idx1,~] = find(SNR(u,b)>=SNRth(:),1,'last');
                CQI(u,b) = single(idx1-1);
                MbpsperRB(u,b) = single(MCS(idx1));
            end
        end
        %% UE-BS association (following the best SNR) and resource allocation
        associatedBSs = zeros(size(SNR,1),NConn);
        [~,associatedBSsaux] = sort(SNR,2,'descend');
        associatedBSs(:,1) = associatedBSsaux(:,1);
        % MC UE list with SNR below the SNR threshold
        for u = 1:Nue
            if SNR(u,associatedBSs(u,1))<SNRmcTh
                associatedBSs(u,2:NConn) = associatedBSsaux(u,2:NConn);
            end
        end
        [MCueList,~]=find(associatedBSs(:,2)>0);
        % Traffic generation for each UE following FTP Model 3
        arrRate = random('Poisson',ArrivalRateAvg,Nue,1);
        requestedRatePerUser = arrRate.*PktSizeAvg;
        for b = 1:Nbs
            [usersAtBS,~] = find(associatedBSs == b);
            usersAtBSRowCol = cell(1);
            for conn = 1:NConn
                [aux,~] = find(associatedBSs(:,conn) == b);
                if sum(CQI(aux,b) == 0)~=0
                    aux(CQI(aux,b)==0) = [];
                end
                usersAtBSRowCol{:,conn} = aux;
            end
            usersAtBS(CQI(usersAtBS,b) == 0) = [];
            if ~isempty(usersAtBS)
                RBperUser = zeros(length(usersAtBS),1);
                % RB requested by UEs depending on the traffic generated and the channel conditions
                ReqRBperUser = ceil(requestedRatePerUser(usersAtBS)./MbpsperRB(usersAtBS,b));
                %% Actual BLER for each UE
                actualBLERperUser_aux = zeros(length(usersAtBS),1);
                for j = 1:length(usersAtBS)
                    auxBler = BLcurves{CQI(usersAtBS(j),b)};
                    auxSNR = find(auxBler(:,1)>=SNR(usersAtBS(j),b),1);
                    if isempty(auxSNR)
                        actualBLERperUser_aux(j) = min(auxBler(:,2)); % Out of BLER curves
                    else
                        actualBLERperUser_aux(j) = auxBler(auxSNR,2);
                    end
                end
                %% Resource scheduling for each UE
                % Weights of RB per UE (i.e., fraction of K)
                weightPerUser = zeros(length(usersAtBS),1);
                switch scheduling
                    case 'RR'
                        weightPerUser = ones(length(usersAtBS),1) / length(usersAtBS);
                    case 'BestCQI'
                        weightPerUser = CQI(usersAtBS,b)./sum(CQI(usersAtBS,b));
                end
                if round(sum(weightPerUser),3)~=1
                    error('The weights sum is not 1');
                end
                auxRBassignment = floor(weightPerUser.*K);
                for rr = 1:length(usersAtBS)
                    if auxRBassignment(rr)>=ReqRBperUser(rr)
                        RBperUser(rr) = ReqRBperUser(rr);
                    else
                        RBperUser(rr) = auxRBassignment(rr);
                    end
                end
                residualRBs = K-sum(RBperUser);
                deltaRBs = ReqRBperUser-RBperUser;
                % This cycle assigns the residual RBs in case the requested RBs per UE are not satisfied
                while residualRBs>0 && any(ReqRBperUser~=RBperUser)
                    [idx,~] = find((ReqRBperUser~=RBperUser)==1);
                    [~,idxSort] = sort(weightPerUser(idx),'descend');
                    weightAux = zeros(length(usersAtBS(idx)),1);
                    switch scheduling
                        case 'RR'
                            weightAux = ones(length(usersAtBS(idx)),1) / length(usersAtBS(idx));
                        case 'BestCQI'
                            weightAux = CQI(usersAtBS(idx),b)./sum(CQI(usersAtBS(idx),b));
                    end
                    if sum(deltaRBs)<=residualRBs
                        RBperUser(idx) = ReqRBperUser(idx);
                    else
                        if sum(floor(weightAux.*residualRBs))==0
                            k = 1;
                            while residualRBs>0
                                if RBperUser(idx(idxSort(k)))<ReqRBperUser(idx(idxSort(k)))
                                    RBperUser(idx(idxSort(k))) = RBperUser(idx(idxSort(k)))+1;
                                    residualRBs = residualRBs-1;
                                end
                                k = k+1;
                                if k>length(idxSort)
                                    k = 1;
                                end
                            end
                        else
                            auxRBres = floor(weightAux.*residualRBs);
                            for rr = 1:length(idx)
                                if auxRBres(rr)<=deltaRBs(idx(rr))
                                    RBperUser(idx(rr)) = RBperUser(idx(rr))+auxRBres(rr);
                                else
                                    RBperUser(idx(rr)) = ReqRBperUser(idx(rr));
                                end
                            end
                        end
                    end
                    residualRBs = K-sum(RBperUser);
                    deltaRBs = ReqRBperUser-RBperUser;
                    if sum(RBperUser)>K
                        error('Assigned more RB than available.');
                    end
                end
                %% Throughput and BLER computation
                listRB(b,:) = [sum(ReqRBperUser) sum(RBperUser)];
                availableRatePerUser_aux = RBperUser .* MbpsperRB(usersAtBS,b);
                throughput_aux = availableRatePerUser_aux .* (1 - actualBLERperUser_aux);
                for conn = 1:NConn
                    actualBLERperUser(usersAtBSRowCol{conn},conn) = actualBLERperUser_aux(1:numel(usersAtBSRowCol{conn}));
                    availableRatePerUser(usersAtBSRowCol{conn},conn) = availableRatePerUser_aux(1:numel(usersAtBSRowCol{conn}));
                    throughput(usersAtBSRowCol{conn},conn) = throughput_aux(1:numel(usersAtBSRowCol{conn}));
                    throughput_aux(1:numel(usersAtBSRowCol{conn})) = [];
                end
            end
        end
        % MC throughput and BLER computation
        for i = 1:Users
            if sum(i==MCueList)
                [~,listConn1] = sort(availableRatePerUser(i,:),'ascend');
                [~,listConn2] = sort(availableRatePerUser(i,:),'descend');
                aux_thr = 0;
                for nconn = 1:NConn
                    if nconn==1
                        aux_thr = throughput(i,listConn1(nconn));
                    else
                        aux_thr = aux_thr.*actualBLERperUser(i,listConn1(nconn)) + throughput(i,listConn1(nconn));
                    end
                end
                throughputMC(i) = aux_thr;
                BLERmc(i) = actualBLERperUser(i,listConn2(1))*actualBLERperUser(i,listConn2(2));
            else
                throughputMC(i) = throughput(i,1);
                BLERmc(i) = actualBLERperUser(i,1);
            end
        end
        % Saving of MC BLERs, BSs usage vs Requested RBs, and MC UE list
        result_throughput{it,l} = single(throughputMC);
        result_bler{it,l} = single(BLERmc);
        result_usage{it,l} = single(listRB);
        result_MCueList{it,l} = single(MCueList);
    end
    if mod(it,100)==0
        disp(strcat('...',num2str(it),'/',num2str(numIterations)));
    end
end
disp(strcat('End (',string(datetime),')'));
toc;

% save('simResultsMultiConn.mat');
