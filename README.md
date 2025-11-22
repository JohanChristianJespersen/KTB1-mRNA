# Overview

Welcome to KTB1-mRNA. This repository contains Simulink files and scripts related to the simulation of a continuous biopharmaceutical process.

# Folders

The project contains three folders, each described in detail below, including the scripts and Simulink files required to run the dynamic model simulations.

- Open-loop simulation | contains the scripts and Simulink file for the open-loop simulation.

- Closed-loop simulation | contains the scripts and Simulink file for the simulation with closed-loop control.

- Closed-loop disturbance simulation | contains the simulation model where disturbances are introduced.

## Files

1. **Simulink**
    - _Final_OpenLoop_ - Open-loop simulation: This Simulink file contains the simulation model for the continuous biopharmaceutical process.
    
    - _Final_Closeloop_ - Closed-loop simulation: This Simulink file contains the simulation model for the continuous biopharmaceutical process with closed-loop control.
      
    - _Final_Closeloop_actuator_ - Actuator & Sensor model: This Simulink file contains actuator and sensor models to increase the realism of the model.
      
    - _Final_Closeloop_five_ - Multiple disturbances: Simulink file introducing Gaussian-distributed variances, step function, ramp signal, a pulse generator that produces a sample-based pulse to the model.
      
    - _Final_sinus_ - Sinusoidal disturbance: This Simulink file contains a sinusoidal disturbance to evaluate the control handle.

2. **Scripts**
    - _CSTR_ - Fermentation Reactor: Script for simulating the continuous reactor.
    - _T_101_ - Buffer tank: Script for simulating Tank T-101 in the process.
    - _T_205_ - Buffer tank: Script for simulating Tank T-205 in the process.
    - _T_302_ - Buffer tank: Script for simulating Tank T-302 in the process.
    - _T_402_ - Buffer tank: Script for simulating Tank T-402 in the process.
    - _T_501_ - Buffer tank: Script for simulating Tank T-501 in the process.
    - _M_501_ - Dynamic mixer: Script for simulating the mixer in the process.
    - _AC_Chromatography_ - Affinity Chromatography: Script for simulating intermediate purification (AC-401:404) in the process.
    - _AEX_Chromatography_ - Anion Exchange Chromatography: Script for simulating initial purification (AEX-201:204) in the process.
    - _AEX_Chromatography2_ - Anion Exchange Chromatography: Script for simulation polishing (AEX-301:304) in the process.
    - _HIC_Chromatography_ - Hydrophobic Interaction Chromatography: Script for simulating polishing (HIC-401:404) based on a two-component system in the process.
    - _Alkalinelysis_ - Alkaline cell lysis: Script for simulation of the alkaline lysis process in a plug flow reactor.
    - _Linearization_ - Linearization: Script for simulation of the pDNA linearization in Reactor R-301.
    - _IVT_ - mRNA transcription: Script for simulation of mRNA transcription in Reactor R-401.
    - _Degradation_reactor_ - DNase I degradation: Script for simulation of DNase I degradation of DNA in Reactor R-402. 
 

# Installation and Usage
**1.** Run the necessary scripts before running the Simulink files.

**2.** Once the scripts are run, you can proceed to run the following Simulink files.

# Additional Information

For any questions or issues please contact jespersen.johan@gmail.com/




Created by Johan Christian Jespersen & Thomas Støvring Sørensen
