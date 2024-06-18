# Performance Evaluation Simulator for Cellular Multi-Connectivity in Emergency Scenarios

This repository contains the MATLAB simulator and necessary files used for the performance evaluation described in the paper titled "_Performance Evaluation of Cellular Multi-Connectivity for Reliable Communications in Emergency Scenarios_" submitted to the **2024 International Conference on Information and Communication Technologies for Disaster Management (ICT-DM 2024)**.

## Overview

The simulator performs a Monte Carlo simulation to compare the performance in terms of BLock Error Rate (BLER) and throughput of single connectivity versus multi-connectivity in a simulated emergency scenario. The evaluation metrics focus on the reliability and robustness of the communication systems under various conditions. Please refer to the paper for the system model description.

## Features

- **Single Connectivity Simulation:** Evaluate performance metrics for traditional single connectivity setups.
- **Multi-Connectivity Simulation:** Assess the advantages of multi-connectivity in terms of reliability and robustness.
- **Monte Carlo Simulation:** Perform extensive simulations to obtain statistically significant results.
- **BLER curves:** The simulator exploits BLER curves as function of SNR obtained from 5G MCS.

## Example

You can perform different analysis changing parameters such as:

- Scenario length
- Resource scheduling technique
- Number of Users
- SNR threshold for MC
- Number of MC links
- BS transmitted power
- Target BLER
- Average packet arrival rate
- Average packet size

## Contributing & Support

Contributions are welcome! Please submit a pull request or open an issue to discuss any changes or improvements. If you want to support our work you can cite the paper, while if you are interested in collaborate with us, please write to **alex.piccioni@univaq.it** or contact the authors of the paper.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This work was supported by University of L'Aquila. The simulator is based on research conducted for the paper submitted to ICT-DM 2024.
