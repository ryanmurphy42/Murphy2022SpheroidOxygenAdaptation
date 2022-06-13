# 2022--Murphy et al.--Growth and adaptation mechanisms of tumour spheroids with time-dependent oxygen availability

Preprint on bioRxiv: https://www.biorxiv.org/content/10.1101/2022.04.24.489294

Experimental data and key MATLAB code used to generate figures to be uploaded June 2022.

Please contact Ryan Murphy for any queries or questions.


Software requirements: 
MATLAB R2021b (v9.11) with the Image Processing Toolbox (v11.4), Optimization Toolbox (v9.2), Global Optimization Toolbox (v4.6), and the Statistics and Machine Learning Toolbox (v12.2). Profile likelihood analysis is performed using code available on the GitHub repository: https://github.com/ryanmurphy42/4DSpheroids_Murphy2021. Parameter estimation, by computing MCMC chains, uses the  \textit{MCMCstat} package available on the GitHub repository https://mjlaine.github.io/mcmcstat/ .

## Guide to using the code and reproducing results

1. Download this GitHub repository
2. Download the MCMCstat package available on the GitHub repository https://mjlaine.github.io/mcmcstat/ .
3. Download the profile likelihood analysis code available on the GitHub repository: https://github.com/ryanmurphy42/4DSpheroids_Murphy2021. 
4. Set 2Code as current folder in MATLAB
5. Open and run each script in 2Code/1FirstAnalysisExperimentalData/ to generate figures of time-evolution of radial measurements and figures for analysis of spheroid snapshots independently. Figures will be saved in the relevant folders in 4Figures and .mat files in 3MATLAB save files.
6. MCMC chains:
	1. Unzip the MCMCstat folder to a folder mcmcstat-master in the same directory as 2Code
	2. Open 2Code/OModels and copy and paste all files from this folder to the mcmcstat-master folder.
	3. Set 2Code as current folder in MATLAB
	4. Open and run each script in 2Code/3MCMC/ to generate MCMC chains. Each script computes in approximately 7 hours on a standard laptop, with the exception of files ending alpha or Rd which compute much quicker due to the simpler model. Figures will be saved in the relevant folders in 4Figures and .mat files in 3MATLAB save files.
