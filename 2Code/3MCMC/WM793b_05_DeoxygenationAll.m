%% 2022--Murphy et al--Growth and adaptation mechanisms of tumour spheroids with time-dependent oxygen availability
% DeoxygenationAll

% This script uses the GitHub package developed by Marko Laine (https://mjlaine.github.io/mcmcstat/)

%% CODE STRUCTURE

% 0) Filepath to save results
% 1) Load data
% 2) Estimate initial outer radius, Ro(0)
% 6) Estimate all parameters
% 6.1) Global minimisation for first guess
% 6.2) Compute mse of model simulated with global min
% 6.3) Run MCMC chains
% 6.4) MCMC Rhat, mean and quantiles of chain
% 6.5) Sample MCMC chain to plot prediction intervals
% 7) Save data to .mat file
% 8) Plots
% 8.0) Load data set from 7
% 8.1) Plot fit to initial outer radius
% 8.4) plot the experimental data
% 8.5) plot the model simulated with the global min vs the experimental data
% 8.6) Plot chains
% 8.7) Plot posterior densities
% 8.9) Plot the mcmc chain error standard deviation
% 8.10) Plot the experimental data v the model simulated with the mean of the chain
% 8.11) Plot the prediction intervals vs experimental data (scatter plot)
% 8.12) Plot the prediction intervals vs experimental data (box chart)
% 9) Change directory back to 2Code

include_plots = 0;

%% 0) Filepath to save results

% open data from 1Data folder
mydir  = pwd; idcs   = strfind(mydir,'\'); newdir = mydir(1:idcs(end)-1);

% create folders
cd ../
if ~exist('3MATLABsavefiles\3MCMC\WM793b\5DeoxygenationAll', 'dir')
    mkdir(['3MATLABsavefiles\3MCMC\WM793b\5DeoxygenationAll']);
end
if ~exist('4Figures\3MCMC\WM793b\5DeoxygenationAll', 'dir')
    mkdir(['4Figures\3MCMC\WM793b\5DeoxygenationAll']);
end


filepath_save = [newdir '\3MATLABsavefiles\3MCMC\WM793b\5DeoxygenationAll\'];
filepath_save_figs =  [newdir '\4Figures\3MCMC\WM793b\5DeoxygenationAll\'];

simulation_id = 'WM793b_05_DeoxygenationAll';

% open directory with the MCMC functions
cd 'mcmcstat-master'

%% 1) Load data

% open data from 1Data folder
filepath_load = [newdir '\3MATLABsavefiles\1FirstAnalysisExperimentalData\WM793b\'];

% remove outliers
WM793b_outer_tmp = load([filepath_load 'WM793b_confocal_scenario_3_outer.mat']);
WM793b_necrotic_tmp = load([filepath_load 'WM793b_confocal_scenario_3_necrotic.mat']);
WM793b_inhibited_tmp = load([filepath_load 'WM793b_confocal_scenario_3_inhibited.mat']);

WM793b_outer.Day = WM793b_outer_tmp.Day(WM793b_outer_tmp.Outlier == 0);
WM793b_outer.Radius = WM793b_outer_tmp.Radius(WM793b_outer_tmp.Outlier == 0);
WM793b_outer.Outlier = WM793b_outer_tmp.Outlier(WM793b_outer_tmp.Outlier == 0);

WM793b_necrotic.Day = WM793b_necrotic_tmp.Day(WM793b_necrotic_tmp.Outlier == 0);
WM793b_necrotic.Radius = WM793b_necrotic_tmp.Radius(WM793b_necrotic_tmp.Outlier == 0);
WM793b_necrotic.Outlier = WM793b_necrotic_tmp.Outlier(WM793b_necrotic_tmp.Outlier == 0);

WM793b_inhibited.Day = WM793b_inhibited_tmp.Day(WM793b_inhibited_tmp.Outlier == 0);
WM793b_inhibited.Radius = WM793b_inhibited_tmp.Radius(WM793b_inhibited_tmp.Outlier == 0);
WM793b_inhibited.Outlier = WM793b_inhibited_tmp.Outlier(WM793b_inhibited_tmp.Outlier == 0);

% load data from normoxia, hypoxia, and deoxygenation alpha and Rd .mat
% files to inform parameter estimation bounds

filepath_load_normoxia = [newdir '\3MATLABsavefiles\3MCMC\WM793b\1Normoxia\'];
filepath_load_hypoxia = [newdir '\3MATLABsavefiles\3MCMC\WM793b\2Hypoxia\'];
filepath_load_alpha = [newdir '\3MATLABsavefiles\3MCMC\WM793b\3DeoxygenationAlpha\'];
filepath_load_Rd  = [newdir '\3MATLABsavefiles\3MCMC\WM793b\4DeoxygenationRd\'];

load_normoxia = load([filepath_load_normoxia 'MCMC_WM793b_01_NormoxiaGreenspan.mat'],"-mat",'params_withnormalfit');
load_hypoxia = load([filepath_load_hypoxia 'MCMC_WM793b_02_HypoxiaGreenspan.mat'],"-mat",'params_withnormalfit');
load_alpha = load([filepath_load_alpha 'WM793b_03_DeoxygenationAlpha.mat'],"-mat",'params_withnormalfit');
load_Rd = load([filepath_load_Rd 'WM793b_04_DeoxygenationRd.mat'],"-mat",'params_withnormalfit');

%% 2) Estimate initial outer radius, Ro(0)

% obtain first day
day_min = min(WM793b_outer.Day);

% select all outer radii data on the first day (and exclude outliers)
OuterRadius_initial = WM793b_outer.Radius(logical((WM793b_outer.Outlier==0).*(WM793b_outer.Day==day_min)));

% fit a normal distribution to data
OuterRadiusInitial_normaldistribution_fit = fitdist(OuterRadius_initial,'Normal');

% output mean and standard deviation
OuterRadiusInitial_normaldistribution_fit.mu;
OuterRadiusInitial_normaldistribution_fit.sigma;


%% 6) Estimate all parameters with MCMC

% theta(1) = alpha_normoxia;
% theta(2) = alpha_hypoxia;
% theta(3) = tau_alpha;
% theta(4) = Rd_normoxia;
% theta(5) = Rd_hypoxia;
% theta(6) = tau_Rd;
% theta(7) = s_normoxia;
% theta(8) = s_hypoxia;
% theta(9) = tau_s;
% theta(10) = lambda_normoxia;
% theta(11) = lambda_hypoxia;
% theta(12) = tau_lambda;
% theta(13) = lambdahat;
% theta(14) = tau_lambdahat;
% theta(15) = outer_radii_at_ts;


data.xdata = [WM793b_outer.Day;WM793b_necrotic.Day;WM793b_inhibited.Day]; % Days
data.ydata = [WM793b_outer.Radius;WM793b_necrotic.Radius;WM793b_inhibited.Radius]; % Outer, Necrotic, Inhibited
variable_id_vec = [ones(length(WM793b_outer.Day),1);2*ones(length(WM793b_necrotic.Day),1);3*ones(length(WM793b_inhibited.Day),1) ];


Omega_unconverted = 3.0318*10^7; %mmHg kg m^(-3)
Omega = Omega_unconverted*(21/160); % % kg m^(-3)

p_infinity = 2; % [percent]

k_unconverted = 2*10^(-9); % [m^2 s^-1]
k = k_unconverted; %*(1/(1*10^(-6))^2)*(1./(60*60*24)); % [micrometre^2 day^(-1)]
t_s=2;

modelfun = @(x,theta) function_deoxygenation_1simulation(theta,...
    k,...
    p_infinity,...
    Omega,...
    t_s,...
    x,...
    max(x),...
    variable_id_vec);

ssfun    = @(theta,data) sum((data.ydata-modelfun(data.xdata,theta)).^2);


%% 6.1) Global minimisation for first guess

% estimate alpha from estimates of Rc from normoxia and hypoxia
p_infinity_21 = 21; % [percent]
p_infinity_2 = 2; % [percent]

alpha_normoxia_first_guess_from_Rc = (6*k*p_infinity_21)./(Omega.*(load_normoxia.params_withnormalfit{4}{2}*1e-6).^2)*1e7; 
alpha_normoxia_lowerbound_from_Rc = (6*k*p_infinity_21)./(Omega.*(load_normoxia.params_withnormalfit{4}{4}*1e-6).^2)*1e7; 
alpha_normoxia_upperbound_from_Rc = (6*k*p_infinity_21)./(Omega.*(load_normoxia.params_withnormalfit{4}{3}*1e-6).^2)*1e7; 

alpha_hypoxia_first_guess_from_Rc = (6*k*p_infinity_2)./(Omega.*(load_hypoxia.params_withnormalfit{4}{2}*1e-6).^2)*1e7; 
alpha_hypoxia_upperbound_from_Rc = (6*k*p_infinity_2)./(Omega.*(load_hypoxia.params_withnormalfit{4}{3}*1e-6).^2)*1e7; 
alpha_hypoxia_lowerbound_from_Rc = (6*k*p_infinity_2)./(Omega.*(load_hypoxia.params_withnormalfit{4}{4}*1e-6).^2)*1e7; 


theta_first_guess = [alpha_normoxia_first_guess_from_Rc,... %alpha_normoxia
    alpha_hypoxia_first_guess_from_Rc,... %alpha_hypoxia
    load_alpha.params_withnormalfit{3}{7},... %tau_alpha
    load_Rd.params_withnormalfit{1}{7},... %Rd_normoxia
    load_Rd.params_withnormalfit{2}{7},... %Rd_hypoxia
    load_Rd.params_withnormalfit{3}{7},... %tau_Rd
    load_normoxia.params_withnormalfit{3}{7},... %s_normoxia
    load_hypoxia.params_withnormalfit{3}{7},... %s_hypoxia
    1,... %tau_s
    load_normoxia.params_withnormalfit{2}{7},... %lambda_normoxia
    load_hypoxia.params_withnormalfit{2}{7},... %lambda_hypoxia
    1,... %tau_lambda
    0.01,... %lambdahat
    1,... %tau_lambdahat
    OuterRadiusInitial_normaldistribution_fit.mu]; %outer_radii_at_ts

fun_likelihood = @(p) ssfun(p,data);
p_first_guess = theta_first_guess;
p_lower_bounds = [max(0.01,alpha_normoxia_lowerbound_from_Rc),...% theta(1) alpha_normoxia
    max(0.01,min(50,alpha_hypoxia_lowerbound_from_Rc)),...% theta(2) = %alpha_hypoxia
    max(0.01,load_alpha.params_withnormalfit{3}{7}-5*load_alpha.params_withnormalfit{3}{8}),...% theta(3) = %tau_alpha
    max(100,min(load_normoxia.params_withnormalfit{1}{7}-5*load_normoxia.params_withnormalfit{1}{8},load_Rd.params_withnormalfit{1}{7}-5*load_Rd.params_withnormalfit{1}{8})),...% theta(4) = %Rd_normoxia
    max(50,min(load_hypoxia.params_withnormalfit{1}{7}-5*load_normoxia.params_withnormalfit{1}{8},load_Rd.params_withnormalfit{2}{7}-5*load_Rd.params_withnormalfit{2}{8})),...% theta(5) = %Rd_hypoxia
    max(0.01,load_Rd.params_withnormalfit{3}{7}-5*load_Rd.params_withnormalfit{3}{8}),...% theta(6) = %tau_Rd
    max(0.01,load_normoxia.params_withnormalfit{3}{7}-5*load_normoxia.params_withnormalfit{3}{8}),...% theta(7) = %s_normoxia
    max(0.01,load_hypoxia.params_withnormalfit{3}{7}-5*load_hypoxia.params_withnormalfit{3}{8}),...% theta(8) = %s_hypoxia
    0.001,...% theta(9) = %tau_s
    max(0.01,load_normoxia.params_withnormalfit{2}{7}-5*load_normoxia.params_withnormalfit{2}{8}),...% theta(10) = %lambda_normoxia
    max(0.01,load_hypoxia.params_withnormalfit{2}{7}-5*load_hypoxia.params_withnormalfit{2}{8}),...% theta(11) = %lambda_hypoxia
    0.001,...% theta(12) = %tau_lambda
    0.001,...% theta(13) = %lambdahat
    0.1,...% theta(14) = %tau_lambdahat
    OuterRadiusInitial_normaldistribution_fit.mu-3*OuterRadiusInitial_normaldistribution_fit.sigma];% theta(15) = %outer_radii_at_ts
p_upper_bounds = [min(25,alpha_normoxia_upperbound_from_Rc),...% theta(1) alpha_normoxia
    min(25,alpha_hypoxia_upperbound_from_Rc),...% theta(2) = %alpha_hypoxia
    load_alpha.params_withnormalfit{3}{7}+5*load_alpha.params_withnormalfit{3}{8},...% theta(3) = %tau_alpha
    max(load_normoxia.params_withnormalfit{1}{7}+5*load_normoxia.params_withnormalfit{1}{8},load_Rd.params_withnormalfit{1}{7}+5*load_Rd.params_withnormalfit{1}{8}),...% theta(4) = %Rd_normoxia
    max(load_hypoxia.params_withnormalfit{1}{7}+5*load_hypoxia.params_withnormalfit{1}{8},load_Rd.params_withnormalfit{2}{7}+5*load_Rd.params_withnormalfit{2}{8}),...% theta(5) = %Rd_hypoxia
    load_Rd.params_withnormalfit{3}{7}+5*load_Rd.params_withnormalfit{3}{8},...% theta(6) = %tau_Rd
    load_normoxia.params_withnormalfit{3}{7}+5*load_normoxia.params_withnormalfit{3}{8},...% theta(7) = %s_normoxia
    load_hypoxia.params_withnormalfit{3}{7}+5*load_hypoxia.params_withnormalfit{3}{8},...% theta(8) = %s_hypoxia
    10,...% theta(9) = %tau_s
    load_normoxia.params_withnormalfit{2}{7}+5*load_normoxia.params_withnormalfit{2}{8},...% theta(10) = %lambda_normoxia
    load_hypoxia.params_withnormalfit{2}{7}+5*load_hypoxia.params_withnormalfit{2}{8},...% theta(11) = %lambda_hypoxia
    10,...% theta(12) = %tau_lambda
    10,...% theta(13) = %lambdahat
    10,...% theta(14) = %tau_lambdahat
    OuterRadiusInitial_normaldistribution_fit.mu+3*OuterRadiusInitial_normaldistribution_fit.sigma];% theta(15) = %outer_radii_at_ts

max_iterations = 100000;
MaxFunctionEvaluations = 100000;
problem = createOptimProblem('fmincon',...
    'objective', fun_likelihood,...
    'x0',p_first_guess,'lb',p_lower_bounds,'ub', p_upper_bounds, 'options',...
    optimoptions(@fmincon,'Algorithm','sqp','Display','off','MaxIterations',max_iterations,'MaxFunctionEvaluations',MaxFunctionEvaluations));

gs = GlobalSearch('Display','iter','NumTrialPoints',5000,'MaxTime',60*60);

tic
rng(14,'twister') % for reproducibility
[p_mle,nLLmin,exitflag_gs,output_gs] = run(gs,problem);
toc

%% 6.2) Compute mse of the model simulated with the global minimum

x_data_outer_pos = 1:length(WM793b_outer.Day);
x_data_necrotic_pos = (length(WM793b_outer.Day)+1):(length(WM793b_outer.Day)+length(WM793b_necrotic.Day));
x_data_inhibited_pos = (length(WM793b_outer.Day)+length(WM793b_necrotic.Day)+1):(length(WM793b_outer.Day)+length(WM793b_necrotic.Day)+length(WM793b_inhibited.Day));

xdata_mse = [data.xdata(x_data_outer_pos);data.xdata(x_data_necrotic_pos);data.xdata(x_data_inhibited_pos)];
variable_id_vec_mse = [ones(length(data.xdata(x_data_outer_pos)),1);2*ones(length(data.xdata(x_data_necrotic_pos)),1);3*ones(length(data.xdata(x_data_inhibited_pos)),1) ];

[ydata_mse] = function_deoxygenation_1simulation(p_mle,...
    k,...
    p_infinity,...
    Omega,...
    t_s,...
    [xdata_mse,xdata_mse,xdata_mse],...
    max(xdata_mse),...
    variable_id_vec_mse);

count_for_mse = length(data.xdata);

% MSE = (1/n)(sumi=1..N (yi-Yi)^2) where yi are observations and Yi is prediction
mse = (1./count_for_mse).*sum((data.ydata - ydata_mse).^2);
save([filepath_save 'MCMC_' simulation_id '_globalmin' '.mat'],'-v7.3');

%% 6.3) Run MCMC


params = {
    {'alpha_normoxia',p_mle(1),p_lower_bounds(1), p_upper_bounds(1),load_alpha.params_withnormalfit{1}{7},3*load_alpha.params_withnormalfit{1}{8}}
    {'alpha_hypoxia',p_mle(2), p_lower_bounds(2), p_upper_bounds(2),load_alpha.params_withnormalfit{2}{7},3*load_alpha.params_withnormalfit{2}{8}}
    {'tau_alpha',p_mle(3),p_lower_bounds(3), p_upper_bounds(3),load_alpha.params_withnormalfit{3}{7},3*load_alpha.params_withnormalfit{3}{8}}
    {'Rd_normoxia',p_mle(4),p_lower_bounds(4), p_upper_bounds(4),load_Rd.params_withnormalfit{1}{7},3*load_Rd.params_withnormalfit{1}{8}}
    {'Rd_hypoxia',p_mle(5),p_lower_bounds(5), p_upper_bounds(5),load_Rd.params_withnormalfit{2}{7},3*load_Rd.params_withnormalfit{2}{8}}
    {'tau_Rd',p_mle(6),p_lower_bounds(6), p_upper_bounds(6),load_Rd.params_withnormalfit{3}{7},3*load_Rd.params_withnormalfit{3}{8}}
    {'s_normoxia',p_mle(7),p_lower_bounds(7), p_upper_bounds(7),load_normoxia.params_withnormalfit{3}{7},3*load_normoxia.params_withnormalfit{3}{8}}
    {'s_hypoxia',p_mle(8),p_lower_bounds(8), p_upper_bounds(8),load_hypoxia.params_withnormalfit{3}{7},3*load_hypoxia.params_withnormalfit{3}{8}}
    {'tau_s',p_mle(9),p_lower_bounds(9), p_upper_bounds(9)}
    {'lambda_normoxia',p_mle(10),p_lower_bounds(10), p_upper_bounds(10),load_normoxia.params_withnormalfit{2}{7},3*load_normoxia.params_withnormalfit{2}{8}}
    {'lambda_hypoxia',p_mle(11),p_lower_bounds(11), p_upper_bounds(11),load_hypoxia.params_withnormalfit{2}{7},3*load_hypoxia.params_withnormalfit{2}{8}}
    {'tau_lambda',p_mle(12),p_lower_bounds(12), p_upper_bounds(12)}
    {'lambdahat',p_mle(13),p_lower_bounds(13), p_upper_bounds(13)}
    {'tau_lambdahat',p_mle(14),p_lower_bounds(14), p_upper_bounds(14)}
    {'outerradii_at_ts',p_mle(15),p_lower_bounds(15), p_upper_bounds(15),p_upper_bounds(15),OuterRadiusInitial_normaldistribution_fit.mu,2*OuterRadiusInitial_normaldistribution_fit.sigma}
    };



model.ssfun  = ssfun;
model.sigma2 = mse; % (initial) error variance from residuals of the global min fit
model.N = length(data.ydata);  % total number of observations
model.S20 = model.sigma2;      % prior mean for sigma2

options.nsimu = 250000;
options.updatesigma = 1;
chain_burnin = 50001;

% run first chain
tic
rng(14,'twister') % for reproducibility
[res_1,chain_tmp_1,s2chain_tmp_1] = mcmcrun(model,data,params,options);
toc
save([filepath_save 'MCMC_' simulation_id '_chain1' '.mat'],'-v7.3');
% run second chain
tic
rng(15,'twister') % for reproducibility
[res_2,chain_tmp_2,s2chain_tmp_2] = mcmcrun(model,data,params,options);
toc
save([filepath_save 'MCMC_' simulation_id '_chain2' '.mat'],'-v7.3');
% run third chain
tic
rng(16,'twister') % for reproducibility
[res_3,chain_tmp_3,s2chain_tmp_3] = mcmcrun(model,data,params,options);
toc
save([filepath_save 'MCMC_' simulation_id '_chain3' '.mat'],'-v7.3');
% run fourth chain
tic
rng(17,'twister') % for reproducibility
[res_4,chain_tmp_4,s2chain_tmp_4] = mcmcrun(model,data,params,options);
toc
save([filepath_save 'MCMC_' simulation_id '_chain4' '.mat'],'-v7.3');

% discard the first 50,000 samples for burn-in
chain_tmp_1a= chain_tmp_1(chain_burnin:options.nsimu,:);
s2chain_tmp_1a = s2chain_tmp_1(chain_burnin:options.nsimu,:);
chain_tmp_2a= chain_tmp_2(chain_burnin:options.nsimu,:);
s2chain_tmp_2a = s2chain_tmp_2(chain_burnin:options.nsimu,:);
chain_tmp_3a= chain_tmp_3(chain_burnin:options.nsimu,:);
s2chain_tmp_3a = s2chain_tmp_3(chain_burnin:options.nsimu,:);
chain_tmp_4a= chain_tmp_4(chain_burnin:options.nsimu,:);
s2chain_tmp_4a = s2chain_tmp_4(chain_burnin:options.nsimu,:);

% combine chains
chain = [chain_tmp_1a
    chain_tmp_2a
    chain_tmp_3a
    chain_tmp_4a];

s2chain = [s2chain_tmp_1a
    s2chain_tmp_2a
    s2chain_tmp_3a
    s2chain_tmp_4a];
res=res_1;

%% 6.4) MCMC Rhat, mean and quantiles of chain

% compute MCMC diagnostic Rhat
Rhat = function_psrf_rhat(chain_tmp_1a,...
    chain_tmp_2a,...
    chain_tmp_3a,...
    chain_tmp_4a);

% compute chain properties
chain_mean = mean(chain);
chain_std = std(chain);
chain_quantiles = quantile(chain,[0.25,0.50,0.75]);

%% 6.5) Sample MCMC chain to plot prediction intervals

time_prediction_interval_one_variable = linspace(2,10)';
len_time_prediction_interval_one_variable = length(time_prediction_interval_one_variable);
time_prediction_interval = [time_prediction_interval_one_variable;time_prediction_interval_one_variable;time_prediction_interval_one_variable];

x_data_outer_pos_prediction_interval = 1:length(time_prediction_interval_one_variable);
x_data_necrotic_pos_prediction_interval = (length(time_prediction_interval_one_variable)+1):(length(time_prediction_interval_one_variable)+length(time_prediction_interval_one_variable));
x_data_inhibited_pos_prediction_interval = (length(time_prediction_interval_one_variable)+length(time_prediction_interval_one_variable)+1):(length(time_prediction_interval_one_variable)+length(time_prediction_interval_one_variable)+length(time_prediction_interval_one_variable));

variable_id_vec_prediction_interval = [ones(len_time_prediction_interval_one_variable,1);2*ones(len_time_prediction_interval_one_variable,1);3*ones(len_time_prediction_interval_one_variable,1) ];

modelfun_prediction_interval = @(x,theta) function_deoxygenation_1simulation(theta,...
    k,...
    p_infinity,...
    Omega,...
    t_s,...
    x,...
    max(x),...
    variable_id_vec_prediction_interval);


% randomly sample from the chain
nsample = 50000;
nsimu  = size(chain,1);

if nsample ==  size(chain,1)
    isample = 1:size(chain,1); % sample whole chain
else
    % random sample from the chain
    isample = ceil(rand(nsample,1)*nsimu);
end

parind = res.parind;
theta  = res.theta;
local  = res.local;

for iisample = 1:nsample
    iisample
    theta(parind) = chain(isample(iisample),:)';
    th  = theta(local==0|local==1);
    y_tmp   = feval(modelfun_prediction_interval,time_prediction_interval,th);
    y = y_tmp + sqrt(s2chain(isample(iisample))).*randn(length(y_tmp),1); % add noise to each value of y
    y(y<0)=0; % enforce non-negative values
    % save
    if iisample == 1
        ysave = zeros(nsample,size(y,1),size(y,2));
    end
    ysave(iisample,:,:) = y;
end



%% 7) Save outputs

close all
save([filepath_save 'MCMC_WM793b_05_DeoxygenationAll' '.mat'],'-v7.3');


%% 8) PLOTS

if include_plots == 1

    %% 8.0) Load experimental data set

    load([filepath_save 'MCMC_WM793b_05_DeoxygenationAll' '.mat']);


    %% 8.1) Plot fit to initial outer radius

    figure(001); clf;
    hold on
    OuterRadius_initial_plot = linspace(OuterRadiusInitial_normaldistribution_fit.mu-3.*OuterRadiusInitial_normaldistribution_fit.sigma,...
        OuterRadiusInitial_normaldistribution_fit.mu+3.*OuterRadiusInitial_normaldistribution_fit.sigma);
    y_OuterRadiusInitial_normaldistribution_fit = pdf(OuterRadiusInitial_normaldistribution_fit,OuterRadius_initial_plot);
    plot(OuterRadius_initial_plot,y_OuterRadiusInitial_normaldistribution_fit,'r')
    xlim([0,400])
    xlabel('Roinit')
    ylabel(['probability'])
    xticks([0,100,200,300,400])
    box on
    grid on
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig001'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig001'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig001'  '.png'])
    print([filepath_save_figs simulation_id '_fig001'],'-depsc2','-painters')






    %% 8.4) Plot experimental data

    figure(004); clf
    hold on
    scatter(data.xdata(x_data_outer_pos),data.ydata(x_data_outer_pos),'g')
    scatter(data.xdata(x_data_necrotic_pos),data.ydata(x_data_necrotic_pos),'k')
    scatter(data.xdata(x_data_inhibited_pos),data.ydata(x_data_inhibited_pos),'m')
    xlim([0,10])
    ylim([0,400])
    xlabel('time [days]')
    ylabel(['Radius'])
    box on
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig004'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig004'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig004'  '.png'])
    print([filepath_save_figs simulation_id '_fig004'],'-depsc2','-painters')

    %% 8.5) plot the model simulated with the global min vs the experimental data

    xdata_globalmin = 2:0.1:10;
    variable_id_vec_globalmin = [ones(length(xdata_globalmin),1);2*ones(length(xdata_globalmin),1);3*ones(length(xdata_globalmin),1) ];
    x_data_outer_pos_globalmin = 1:length(xdata_globalmin);
    x_data_necrotic_pos_globalmin = (length(xdata_globalmin)+1):(length(xdata_globalmin)+length(xdata_globalmin));
    x_data_inhibited_pos_globalmin = (length(xdata_globalmin)+length(xdata_globalmin)+1):(length(xdata_globalmin)+length(xdata_globalmin)+length(xdata_globalmin));

    [ydata_globalmin] = function_deoxygenation_1simulation(p_mle,...
        k,...
        p_infinity,...
        Omega,...
        t_s,...
        [xdata_globalmin,xdata_globalmin,xdata_globalmin],...
        max(xdata_globalmin),...
        variable_id_vec_globalmin);

    figure(005); clf
    hold on
    scatter(data.xdata(x_data_outer_pos),data.ydata(x_data_outer_pos),'g')
    scatter(data.xdata(x_data_necrotic_pos),data.ydata(x_data_necrotic_pos),'k')
    scatter(data.xdata(x_data_inhibited_pos),data.ydata(x_data_inhibited_pos),'m')
    hold on
    plot(xdata_globalmin,ydata_globalmin(x_data_outer_pos_globalmin),'g--')
    plot(xdata_globalmin,ydata_globalmin(x_data_necrotic_pos_globalmin),'k--')
    plot(xdata_globalmin,ydata_globalmin(x_data_inhibited_pos_globalmin),'m--')
    xlim([0,10])
    ylim([0,500])
    xlabel('time [days]')
    ylabel(['Radius'])
    box on

    saveas(gcf, [filepath_save_figs  simulation_id  '_fig005'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig005'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig005'  '.png'])
    print([filepath_save_figs simulation_id '_fig005'],'-depsc2','-painters')




    %% 8.6) Plot the chains

    figure(006); clf
    mcmcplot(chain,[],res,'chainpanel');
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig006'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig006'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig006'  '.png'])
    print([filepath_save_figs simulation_id '_fig006'],'-depsc2','-painters')

    %% 8.7) Plot the posterior densities

    figure(007); clf
    mcmcplot(chain,[],res,'denspanel');
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig007'  '.fig'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig007'  '.pdf'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig007'  '.png'])
    print([filepath_save_figs simulation_id '_fig007'],'-depsc2','-painters')


    figure(008); clf;
    mcmcplot_density(chain,res,filepath_save_figs,simulation_id);
    saveas(gcf, [filepath_save_figs simulation_id  '_fig008' '.fig'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig008' '.pdf'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig008' '.png'])
    print([filepath_save_figs simulation_id '_fig008'],'-depsc2','-painters')

    %% 8.9) Plot the mcmc chain error standard deviation

    figure(010); clf
    mcmcplot(sqrt(s2chain),[],[],'hist')
    title('Error std posterior')

    % add prior distribution to the plot
    if res.N0>0
        xl = xlim; xx = linspace(xl(1),xl(2));
        hold on
        plot(xx,invchi1pf(xx,res.N0,sqrt(res.S20)))
        hold off
        legend('posterior','prior')
    end
    xlim([0,40])
    saveas(gcf, [filepath_save_figs   simulation_id '_fig010'  '.fig'])
    saveas(gcf, [filepath_save_figs   simulation_id '_fig010'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig010'  '.png'])
    print([filepath_save_figs simulation_id '_fig010'],'-depsc2','-painters')


    %% 8.10) Plot the experimental data v the model simulated with the mean of the chain

    xdata_meanchain = 2:0.1:10;
    variable_id_vec_meanchain = [ones(length(xdata_meanchain),1);2*ones(length(xdata_meanchain),1);3*ones(length(xdata_meanchain),1) ];

    [ydata_meanchain] = function_deoxygenation_1simulation(mean(chain),...
        k,...
        p_infinity,...
        Omega,...
        t_s,...
        [xdata_meanchain,xdata_meanchain,xdata_meanchain],...
        max(xdata_meanchain),...
        variable_id_vec_meanchain);

    x_data_outer_pos_meanchain = 1:length(xdata_meanchain);
    x_data_necrotic_pos_meanchain = (length(xdata_meanchain)+1):(length(xdata_meanchain)+length(xdata_meanchain));
    x_data_inhibited_pos_meanchain = (length(xdata_meanchain)+length(xdata_meanchain)+1):(length(xdata_meanchain)+length(xdata_meanchain)+length(xdata_meanchain));

    figure(011); clf
    hold on
    scatter(data.xdata(x_data_outer_pos),data.ydata(x_data_outer_pos),'g')
    scatter(data.xdata(x_data_necrotic_pos),data.ydata(x_data_necrotic_pos),'k')
    scatter(data.xdata(x_data_inhibited_pos),data.ydata(x_data_inhibited_pos),'m')
    hold on
    plot(xdata_meanchain,ydata_meanchain(x_data_outer_pos_meanchain),'g--')
    plot(xdata_meanchain,ydata_meanchain(x_data_necrotic_pos_meanchain),'k--')
    plot(xdata_meanchain,ydata_meanchain(x_data_inhibited_pos_meanchain),'m--')
    xlim([0,10])
    ylim([0,400])
    xlabel('time [days]')
    ylabel(['Radius'])
    saveas(gcf, [filepath_save_figs   simulation_id '_fig011'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig011'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig011'  '.png'])
    print([filepath_save_figs simulation_id '_fig011'],'-depsc2','-painters')

    %%  8.11) Plot the prediction intervals vs experimental data (scatter plot)

    % for each time point, estimate the model values at the quantiles
    lims = [0.005,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.995];
    predlims = plims(ysave(:,:,1),lims);

    % PLOT
    figure(012); clf
    np = size(predlims,1);
    nn = (np+1)/2; % median
    np = nn-1;
    % OUTER
    dimc_outer = [0.8 1.0 0.8]; % dimmest (lightest) color OUTER- green
    fillyy(time_prediction_interval(x_data_outer_pos_prediction_interval),predlims(1,x_data_outer_pos_prediction_interval),predlims(2*nn-1,x_data_outer_pos_prediction_interval),dimc_outer);
    hold on
    for k=2:(nn-1)
        fillyy(time_prediction_interval(x_data_outer_pos_prediction_interval),predlims(k,x_data_outer_pos_prediction_interval),predlims(2*nn-k,x_data_outer_pos_prediction_interval),dimc_outer.*0.9.^(k-1));
    end

    % INHIBITED
    dimc_inhibited = [1.0 0.8 0.9]; % dimmest (lightest) color INHIBITED - magenta
    fillyy(time_prediction_interval(x_data_inhibited_pos_prediction_interval),predlims(1,x_data_inhibited_pos_prediction_interval),predlims(2*nn-1,x_data_inhibited_pos_prediction_interval),dimc_inhibited);
    hold on
    for k=2:(nn-1)
        fillyy(time_prediction_interval(x_data_inhibited_pos_prediction_interval),predlims(k,x_data_inhibited_pos_prediction_interval),predlims(2*nn-k,x_data_inhibited_pos_prediction_interval),dimc_inhibited.*0.9.^(k-1));
    end

    % NECROTIC
    dimc_necrotic = [0.8 0.8 0.8]; % dimmest (lightest) color NECROTIC - grey
    fillyy(time_prediction_interval(x_data_necrotic_pos_prediction_interval),predlims(1,x_data_necrotic_pos_prediction_interval),predlims(2*nn-1,x_data_necrotic_pos_prediction_interval),dimc_necrotic);
    hold on
    for k=2:(nn-1)
        fillyy(time_prediction_interval(x_data_necrotic_pos_prediction_interval),predlims(k,x_data_necrotic_pos_prediction_interval),predlims(2*nn-k,x_data_necrotic_pos_prediction_interval),dimc_necrotic.*0.9.^(k-1));
    end

    % PLOT experimental data
    hold on
    scatter(data.xdata(x_data_outer_pos),data.ydata(x_data_outer_pos),'g')
    scatter(data.xdata(x_data_necrotic_pos),data.ydata(x_data_necrotic_pos),'k')
    scatter(data.xdata(x_data_inhibited_pos),data.ydata(x_data_inhibited_pos),'m')

    % figure properties
    xlim([0,10])
    ylim([0,400])
    box on
    xlabel('time [days]')
    ylabel('Radius [microns]')
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig012'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig012'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig012'  '.png'])
    print([filepath_save_figs simulation_id '_fig012'],'-depsc2','-painters')

    %%  8.12) Plot the prediction intervals vs experimental data (box charts)

    fig_boxchart_predinterval = openfig([newdir '\4Figures\1FirstAnalysisExperimentalData\WM793b\Scenario_3_OIN_1.fig']);

    figure(fig_boxchart_predinterval)
    hold on
    np = size(predlims,1);
    nn = (np+1)/2; % median
    np = nn-1;
    % OUTER
    dimc_outer = [0.8 1.0 0.8]; % dimmest (lightest) color OUTER- green
    fillyy(time_prediction_interval(x_data_outer_pos_prediction_interval),predlims(1,x_data_outer_pos_prediction_interval),predlims(2*nn-1,x_data_outer_pos_prediction_interval),dimc_outer);
    hold on
    for k=2:(nn-1)
        fillyy(time_prediction_interval(x_data_outer_pos_prediction_interval),predlims(k,x_data_outer_pos_prediction_interval),predlims(2*nn-k,x_data_outer_pos_prediction_interval),dimc_outer.*0.9.^(k-1));
    end

    % INHIBITED
    dimc_inhibited = [1.0 0.8 0.9]; % dimmest (lightest) color INHIBITED - magenta
    fillyy(time_prediction_interval(x_data_inhibited_pos_prediction_interval),predlims(1,x_data_inhibited_pos_prediction_interval),predlims(2*nn-1,x_data_inhibited_pos_prediction_interval),dimc_inhibited);
    hold on
    for k=2:(nn-1)
        fillyy(time_prediction_interval(x_data_inhibited_pos_prediction_interval),predlims(k,x_data_inhibited_pos_prediction_interval),predlims(2*nn-k,x_data_inhibited_pos_prediction_interval),dimc_inhibited.*0.9.^(k-1));
    end

    % NECROTIC
    dimc_necrotic = [0.8 0.8 0.8]; % dimmest (lightest) color NECROTIC - grey
    fillyy(time_prediction_interval(x_data_necrotic_pos_prediction_interval),predlims(1,x_data_necrotic_pos_prediction_interval),predlims(2*nn-1,x_data_necrotic_pos_prediction_interval),dimc_necrotic);
    hold on
    for k=2:(nn-1)
        fillyy(time_prediction_interval(x_data_necrotic_pos_prediction_interval),predlims(k,x_data_necrotic_pos_prediction_interval),predlims(2*nn-k,x_data_necrotic_pos_prediction_interval),dimc_necrotic.*0.9.^(k-1));
    end


    % figure properties
    xlim([0,10])
    ylim([0,400])
    box on
    xlabel('time [days]')
    ylabel('Radius [microns]')


    saveas(gcf, [filepath_save_figs  simulation_id  '_fig013_BoxchartOIN_with_predictionintervals'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig013_BoxchartOIN_with_predictionintervals'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig013_BoxchartOIN_with_predictionintervals'  '.png'])
    print([filepath_save_figs simulation_id '_fig013_BoxchartOIN_with_predictionintervals'],'-depsc2','-painters')




end

%% 9) Change directory back to 2Code
cd ../
cd 2Code
