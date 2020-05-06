# Ba,Sr,Ra-SO4_Co-precipitation
###### Simulation of precipitation processes during mixing of two solutions

## Table of Contents:
- [Installation](#installation)
- [Project Motivation](#project-motivation)
- [File Descriptions](#file-descriptions)
- [Results](#results)
- [Licensing, Authors, Acknowledgements](#end)

## Installation
There is no necessary libraries to run the code here beyond the Matlab default settings. The code should run with no issues using Matlab program.

## Project Motivation
With this project, I was interestested in examining the co-precipitation behavior of (Ba,Sr,Ra)SO4 by mixing two solutions:
- Solution 1 with different combinations of Ba, Sr, Ra concentrations, and
- Solution 2 with SO4
in three types of background:
- Freshwater condition with no added salts beyond Cl for Solution 1 and Na for Solution 2
- 2 M NaCl, and
- 2 M NaCl plus 0.5 M CaCl2 to simulate saline and brine environments

## File Descriptions
There are 8 .m matlab files available to showcase work to address the above questions. 

The `Main_program_publish_0509.m` is the main file that runs the program repository and return co-precipitation results.
The `gamma_fresh.m`, `gamma_pitzer.m`, and `gamma_sit.m` function-files provide three methods for computing activity coefficients for involved ions under freshwater (`gamma_fresh.m`) and saline conditions (`gamma_pitzer.m` and `gamma_sit.m`).
The `equ_noBarite.m`, `equ_noCeles.m`, and `equation_sets.m` function-files provide multiple equations to conclude three types of possible co-precipitation procedure happen in one mini-step.
The `Kd_coeff.m` function-file will generate the partition coefficients for `Ra` in `BaSO4` and `SrSO4`, `Ba` in `SrSO4`, and `Sr` in `BaSO4`.

## Results
By choosing the types of Ba-Sr-Ra combination and reaction background, the project will generate a table summarizing changes in ionic concentrations, activity coefficients, partition coefficients, and solid compositions in 10000 steps of titration of 1 L Solution 1 into 100 mL Solution 2.

The main findings of the project can be found at the thesis Chapter 5 available [here](https://search.proquest.com/openview/6ecafcfa25372299a3e4a52b61708830/1?pq-origsite=gscholar&cbl=18750&diss=y).
<a name="end"></a>
## Licensing, Authors, Acknowledgements
Credit is given to the advisors Prof. Xiahong Feng, Prof. Eric Posmentier, Prof. Leslie Sonder, and Prof. Devon Renock of the Earth Sciences Department at Dartmouth College. Otherwise, feel free to use the code here as you would like!

