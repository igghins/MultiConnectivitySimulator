# Performance Evaluation Simulator for Cellular Multi-Connectivity in Emergency Scenarios

This repository contains the MATLAB simulator and necessary files used for the performance evaluation described in the paper titled "_Performance Evaluation of Multi-Connectivity for Massive-URLLC in Emergency Scenarios_" submitted to the **2024 International Conference on Information and Communication Technologies for Disaster Management (ICT-DM 2024)**.

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

Contributions are welcome! Please submit a pull request or open an issue to discuss any changes or improvements. If you want to cite the simulator, please use:

- A. Piccioni, A. Marotta, C. Rinaldi and F. Graziosi, "_Enhancing Mobile Networks for Urban Air Mobility Connectivity,_" in **IEEE Networking Letters**, vol. 6, no. 2, pp. 110-114, June 2024, doi: 10.1109/LNET.2024.3390610.
- A. Piccioni, A. Marotta, P. Di Marco and F. Graziosi, "_Performance Evaluation of Multi-Connectivity for Massive-URLLC in Emergency Scenarios_", ACCEPTED FOR PUBLICATION IN **2024 International Conference on Information and Communication Technologies for Disaster Management (ICT-DM 2024)**.

If you are interested in a collaboration, please write me (**alex.piccioni@univaq.it**) or contact my co-authors.

## License

This project is licensed under the MIT License - see the [MIT](https://choosealicense.com/licenses/mit/) file for details.

## Acknowledgments

This work was supported by University of L'Aquila.
