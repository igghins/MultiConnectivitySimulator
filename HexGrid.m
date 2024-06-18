function [locationMatrix] = HexGrid(origin,ISD,L,BSheight)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    localOrigin = [0 0];
    xMax = localOrigin(1)+L/2-1;
    xMin = localOrigin(1)-L/2;
    yMax = localOrigin(2)+L/2-1;
    yMin = localOrigin(2)-L/2;

    verticalDistance = ISD; % diagonale cerchio interno
    horizontalDistance = sqrt(3)*verticalDistance; % 

    %% horizontal coordinates
    % 1/2 number - 1 of hexagonal structures in the row that contains (0,0)
    nHexagonsRowEvenMax = floor(xMax/horizontalDistance);
    nHexagonsRowEvenMin =  ceil(xMin/horizontalDistance);
    % 1/2 number of hexagonal structures in the rows that are shifted
    nHexagonsRowOddMax = floor((xMax+0.5*horizontalDistance)/horizontalDistance);
    nHexagonsRowOddMin =  ceil((xMin+0.5*horizontalDistance)/horizontalDistance);
    
    % biggest x coordinate for even rows
    xMaxEven = nHexagonsRowEvenMax*horizontalDistance;
    xMinEven = nHexagonsRowEvenMin*horizontalDistance;
    xLocationEven = xMinEven:horizontalDistance:xMaxEven;
    
    % biggest x coordinate for odd rows
    xMaxOdd = (nHexagonsRowOddMax-0.5)*horizontalDistance;
    xMinOdd = (nHexagonsRowOddMin-0.5)*horizontalDistance;
    xLocationOdd = xMinOdd:horizontalDistance:xMaxOdd;
    
    %% vertical coordinates
    % 1/2 number - 1 of hexagonal structures in the column that contains (0,0)
    nHexagonsColumnEvenMax = floor(yMax/verticalDistance);
    nHexagonsColumnEvenMin =  ceil(yMin/verticalDistance);
    % 1/2 number of hexagonal structures in the rows that are 'verschoben'-
    nHexagonsColumnOddMax = floor((yMax+0.5*verticalDistance)/verticalDistance);
    nHexagonsColumnOddMin =  ceil((yMin+0.5*verticalDistance)/verticalDistance);
    
    % biggest y coordinate for even rows
    yMaxEven = nHexagonsColumnEvenMax*verticalDistance;
    yMinEven = nHexagonsColumnEvenMin*verticalDistance;
    yLocationEven = yMinEven:verticalDistance:yMaxEven;
    
    % biggest y coordinate for odd rows
    yMaxOdd = (nHexagonsColumnOddMax-0.5)*verticalDistance;
    yMinOdd = (nHexagonsColumnOddMin-0.5)*verticalDistance;
    yLocationOdd = yMinOdd:verticalDistance:yMaxOdd;
    
    % create grid
    [X1, Y1] = meshgrid(xLocationEven,yLocationEven);
    [X2, Y2] = meshgrid(xLocationOdd,yLocationOdd);
    locationMatrix = [[X1(:); X2(:)] [Y1(:);Y2(:)] BSheight*ones(size([X1(:); X2(:)]))];
    locationMatrix(:,1:2) = locationMatrix(:,1:2)+origin;
    for i = 1:size(locationMatrix,1)
        if locationMatrix(i,1)>L/2
            locationMatrix(i,1) = locationMatrix(i,1)-L;
        end
        if locationMatrix(i,2)>L/2
            locationMatrix(i,2) = locationMatrix(i,2)-L;
        end
    end
end

