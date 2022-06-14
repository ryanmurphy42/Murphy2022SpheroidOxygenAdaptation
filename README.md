# Murphy et al. (2022) Growth and adaptation mechanisms of tumour spheroids with time-dependent oxygen availability

Preprint on bioRxiv: https://www.biorxiv.org/content/10.1101/2022.04.24.489294

This repository holds experimental data and key MATLAB code used to generate figures in the manuscript.

Please contact Ryan Murphy for any queries or questions.

Software requirements: 
MATLAB R2021b (v9.11) with the Image Processing Toolbox (v11.4), Optimization Toolbox (v9.2), Global Optimization Toolbox (v4.6), and the Statistics and Machine Learning Toolbox (v12.2). Profile likelihood analysis is performed using code available on the GitHub repository: https://github.com/ryanmurphy42/4DSpheroids_Murphy2021. Parameter estimation, by computing MCMC chains, uses the MCMCstat package available on the GitHub repository https://mjlaine.github.io/mcmcstat/ .

## Guide to using the code

1. Download this GitHub repository
2. Download the MCMCstat package available on the GitHub repository https://mjlaine.github.io/mcmcstat/ .
3. Download the profile likelihood analysis code available on the GitHub repository: https://github.com/ryanmurphy42/4DSpheroids_Murphy2021. 
4. Set 2Code as current folder in MATLAB
5. Open and run each script in 2Code/1FirstAnalysisExperimentalData/ to generate figures of time-evolution of radial measurements and figures for analysis of spheroid snapshots independently. Figures will be saved in the relevant folders in 4Figures and .mat files in 3MATLABsavefiles.
6. Profile likelihood analysis
	1. Unzip the profile likelihood analysis code available on the GitHub repository: https://github.com/ryanmurphy42/4DSpheroids_Murphy2021 and copy files to 2Code/3Profilelikelihood, and do not replace function_load_simulation_settings.m. Note that function_load_simulation_settings.m contain simulation settings for these experiments.
	2. Copy scenario 1 and scenario 2 .mat files from 3MATLABsavefiles\1FirstAnalysisExperimentalData to 2Code/3Profilelikelihood
	3. Use code as explained on https://github.com/ryanmurphy42/4DSpheroids_Murphy2021 
7. MCMC chains:
	1. Unzip the MCMCstat folder to a folder mcmcstat-master in the same directory as 2Code
	2. Open 2Code/OModels and copy and paste all files from this folder to the mcmcstat-master folder.
	3. Set 2Code as current folder in MATLAB
	4. Open and run each script in 2Code/3MCMC/ to generate MCMC chains. Each script computes in approximately 7 hours on a standard laptop, with the exception of files ending alpha or Rd which compute much quicker due to the simpler model. Figures will be saved in the relevant folders in 4Figures and .mat files in 3MATLAB save files.
