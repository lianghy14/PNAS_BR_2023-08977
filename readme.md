<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>

<!-- ABOUT THE PROJECT -->
## About The Simulation

The Seoul Halloween crowd-crush disaster, which occurred in an alleyway of the Itaewon neighborhood in Yongsan district in Seoul, Korea on October 29, 2022, resulted in the deaths of 159 people, comprising 133 Korean nationals and 26 foreign nationals.

According to the official report, over 300 victims were concentrated in an $18.24$ $\mathrm{m^{2}}$ ($5.7$ $\mathrm{m}$ $\times$ $3.2$ $\mathrm{m}$) area of the alleyway. A conceptual layout of the alleyway and the aforementioned area are depicted as follows.

![alt text](https://github.com/lianghy14/PNAS_BR_2023-08977/blob/39298e0c41d4482d41459d507d36b84901703982/Figures/fig_concept.png)

Based on empirical analysis, the numerical simulation is performed over a $104\times61$ $\mathrm{m^2}$ `$\mathbf{H}$'-shaped area, in which there are six pedestrian streams: two opposing streams in the south street, two in the north street, and two in the alleyway. They were assigned different boundary conditions, i.e., for inflow boundaries $\Gamma_O^{(k)}$, outflow boundaries $\Gamma_D^{(k)}$, and solid boundaries $\Gamma_H^{(k)}$. Six pedestrian streams are considered in the simulation geometry (Layout size: $104\times61$ $\mathrm{m^2}$ mesh size: $208\times122$).

![alt text](https://github.com/lianghy14/PNAS_BR_2023-08977/blob/9030648d4a6f1568b46adebfa908df0fb239ea01/Figures/fig_layout1.png)
![alt text](https://github.com/lianghy14/PNAS_BR_2023-08977/blob/9030648d4a6f1568b46adebfa908df0fb239ea01/Figures/fig_layout2.png)

This project includes all source codings and simulation results of the Seoul Halloween crowd-crush. Since the orginal data is too large, only visualized plots of density at every 60 seconds is included in the `Results`.

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

The codings run on MATLAB with FORTRAN complier. The following software is required to run the simulation.
* Matlab R2022a
* Visual Studio Community 2022
* Intel oneAPI Base Toolkit 2023
* Intel oneAPI HPC Toolkit 2023

<!-- USAGE EXAMPLES -->
## Usage

Run the `MAT_TVDRK2.m` in MATLAB to start the simulation. Simulation results will be recorded in `Results`.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

Howie LIANG - lianghy@connect.hku.hk

Project Link: [https://github.com/lianghy14/PNAS_BR_2023-08977/blob/main](https://github.com/lianghy14/PNAS_BR_2023-08977/blob/main)

<p align="right">(<a href="#readme-top">back to top</a>)</p>
