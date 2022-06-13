%% 2022--Murphy et al--Growth and adaptation mechanisms of tumour spheroids with time-dependent oxygen availability
% ReoxygenationRd

% This script uses the GitHub package developed by Marko Laine (https://mjlaine.github.io/mcmcstat/)

%% CODE STRUCTURE

% 0) Filepath to save results
% 1) Load data
% 1.1) Estimate mean and standard deviation of Rd_normoxia
% 1.2) Estimate mean and standard deviation of Rd_hypoxia
% 2) Estimate all parameters (Rd_normoxia,Rd_hypoxia, tau_Rd)
% 2.1) Global minimisation for first guess
% 2.2) Compute mse of model simulated with global min
% 2.3) Run MCMC chains
% 3) Save data to .mat file
% 4) Plots
% 4.0) Load data set from 7
% 4.1) plot the experimental data
% 4.2) plot the model simulated with the global min vs the experimental data
% 4.3) Plot chains
% 4.4) Plot posterior densities
% 4.5) Plot MCMC pairs
% 4.6) Plot the mcmc chain error standard deviation
% 4.7) Plot the experimental data v the model simulated with the mean of the chain
% 4.8) Plot the prediction intervals vs experimental data (scatter plot)

include_plots = 1;

%% 0) Filepath to save results

% open data from 1Data folder
mydir  = pwd; idcs   = strfind(mydir,'\'); newdir = mydir(1:idcs(end)-1);

% create folders
cd ../
if ~exist('3MATLABsavefiles\3MCMC\WM164\8ReoxygenationRd', 'dir')
    mkdir(['3MATLABsavefiles\3MCMC\WM164\8ReoxygenationRd']);
end
if ~exist('4Figures\3MCMC\WM164\8ReoxygenationRd', 'dir')
    mkdir(['4Figures\3MCMC\WM164\8ReoxygenationRd']);
end

filepath_save = [newdir '\3MATLABsavefiles\3MCMC\WM164\8ReoxygenationRd\'];
filepath_save_figs =  [newdir '\4Figures\3MCMC\WM164\8ReoxygenationRd\'];

simulation_id = 'WM164_08_ReoxygenationRd';

% open directory with the MCMC functions
cd 'mcmcstat-master'

%% 1) Load data

% open data from 1Data folder
filepath_load = [newdir '\3MATLABsavefiles\1FirstAnalysisExperimentalData\WM164\'];

reoxygenation_Rd_data= load([filepath_load 'WM164_oxygenwasteanalysis_reoxy.mat']);

data.xdata = reoxygenation_Rd_data.confocal_WM164.AllData.Day(reoxygenation_Rd_data.Rd_vector>=0); % select all data and filter to valid data (i.e when Rd >=0)
data.ydata = reoxygenation_Rd_data.Rd_vector(reoxygenation_Rd_data.Rd_vector>=0);

%% 1.1) Estimate mean and standard deviation of Rd_normoxia

% obtain first day
day_min = min(data.xdata);

% select all estimates of Rd on the first day (and exclude outliers)
Rd_initial_time = data.ydata(logical(data.xdata==day_min));

% fit a normal distribution to data
Rd_initial_time_normaldistribution_fit = fitdist(Rd_initial_time,'Normal');

% output mean and standard deviation
Rd_initial_time_normaldistribution_fit.mu;
Rd_initial_time_normaldistribution_fit.sigma;


%% 1.2) Estimate mean and standard deviation of Rd_hypoxia

% obtain last day
day_max = max(data.xdata);

% select all estimates of Rd on the first day (and exclude outliers)
Rd_final_time = data.ydata(logical(data.xdata==day_max));

% fit a normal distribution to data
Rd_final_time_normaldistribution_fit = fitdist(Rd_final_time,'Normal');

% output mean and standard deviation
Rd_final_time_normaldistribution_fit.mu;
Rd_final_time_normaldistribution_fit.sigma;

%% 2) Estimate all parameters with MCMC (Rd_normoxia, Rd_hypoxia, tau_Rd)

t_s = 2; % time when switches from normoxia to hypoxia

% theta(1) %Rd_normoxia,...
% theta(2) %Rd_hypoxia,...
% theta(3) %tau_Rd

modelfun = @(x,theta) theta(1) + (theta(2)-theta(1)).*exp(-(1/theta(3))*(x-t_s));
ssfun    = @(theta,data) sum((data.ydata-modelfun(data.xdata,theta)).^2);

%% 2.1) Global minimisation for first guess


theta_first_guess =  [Rd_final_time_normaldistribution_fit.mu; % Rd_normoxia
    Rd_initial_time_normaldistribution_fit.mu;  % Rd_hypoxia
    1]; % tau_Rd

fun_likelihood = @(p) ssfun(p,data);

p_first_guess = theta_first_guess;
p_lower_bounds =  [max(0,Rd_final_time_normaldistribution_fit.mu - 5*Rd_final_time_normaldistribution_fit.sigma); % Rd_normoxia
    max(0, Rd_initial_time_normaldistribution_fit.mu - 5*Rd_initial_time_normaldistribution_fit.sigma);  % Rd_hypoxia
    0]; % tau_Rd
p_upper_bounds =  [Rd_final_time_normaldistribution_fit.mu + 5*Rd_final_time_normaldistribution_fit.sigma; % Rd_normoxia
    Rd_initial_time_normaldistribution_fit.mu + 5*Rd_initial_time_normaldistribution_fit.sigma;  % Rd_hypoxia
    100]; % tau_Rd

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

%% 2.2) Compute mse of the model simulated with the global minimum

xdata_mse = [data.xdata];
[ydata_mse] = modelfun(xdata_mse,p_mle);
count_for_mse = length(data.xdata);

% MSE = (1/n)(sumi=1..N (yi-Yi)^2) where yi are observations and Yi is prediction
mse = (1./count_for_mse).*sum((data.ydata - ydata_mse).^2);


%% 2.3) Run MCMC

params = {
    {'Rd normoxia', p_mle(1), p_lower_bounds(1),p_upper_bounds(1),Rd_final_time_normaldistribution_fit.mu,3*Rd_final_time_normaldistribution_fit.sigma}
    {'Rd hypoxia', p_mle(2), p_lower_bounds(2),p_upper_bounds(2),Rd_initial_time_normaldistribution_fit.mu,3*Rd_initial_time_normaldistribution_fit.sigma}
    {'tau Rd', p_mle(3), p_lower_bounds(3),p_upper_bounds(3)}
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
% run second chain
tic
rng(15,'twister') % for reproducibility
[res_2,chain_tmp_2,s2chain_tmp_2] = mcmcrun(model,data,params,options);
toc
% run third chain
tic
rng(16,'twister') % for reproducibility
[res_3,chain_tmp_3,s2chain_tmp_3] = mcmcrun(model,data,params,options);
toc
% run fourth chain
tic
rng(17,'twister') % for reproducibility
[res_4,chain_tmp_4,s2chain_tmp_4] = mcmcrun(model,data,params,options);
toc

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

% compute MCMC diagnostic Rhat
Rhat = function_psrf_rhat(chain_tmp_1a,...
    chain_tmp_2a,...
    chain_tmp_3a,...
    chain_tmp_4a);

% compute chain properties
chain_mean = mean(chain);
chain_std = std(chain);
chain_quantiles = quantile(chain,[0.25,0.50,0.75]);

% Estimate mean and standard deviation to MCMC chain parameters to inform parameter bounds in deoxygenation, re-oxygenation simulations

params_withnormalfit = params;

for r=1:size(params_withnormalfit,1)
    params_withnormalfit{r}{7} = chain_mean(r);
    params_withnormalfit{r}{8} = chain_std(r);
end



%% 3) Save data to .mat file

save([filepath_save 'WM164_08_ReoxygenationRd' '.mat'],'-v7.3');

%% 4) Plots

if include_plots == 1

    %% 4.1) plot the experimental data

    figure(001); clf
    plot(data.xdata,data.ydata,'s');
    xlim([0,10])
    xlabel('time [days]')
    ylabel('Rd [micrometres]')
    saveas(gcf, [filepath_save_figs simulation_id   '_fig001'  '.fig'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig001'  '.pdf'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig001'  '.png'])
    print([filepath_save_figs  '_fig001'],'-depsc2','-painters')

    %% 4.2) plot the model simulated with the global min vs the experimental data


    xdata_globalmin = 2:0.1:10;
    [ydata_globalmin] = modelfun(xdata_globalmin,p_mle);
    figure(002); clf
    hold on
    plot(data.xdata,data.ydata,'s');
    plot(xdata_globalmin,ydata_globalmin)
    xlim([0,10])
    xlabel('time [days]')
    ylabel('Rd [micrometres]')
    box on

    saveas(gcf, [filepath_save_figs  simulation_id  '_fig002'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig002'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig002'  '.png'])
    print([filepath_save_figs simulation_id '_fig002'],'-depsc2','-painters')


    %%  4.3) Plot chains

    figure(003); clf
    mcmcplot(chain,[],res,'chainpanel');
    saveas(gcf, [filepath_save_figs simulation_id   '_fig003'  '.fig'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig003'  '.pdf'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig003'  '.png'])
    print([filepath_save_figs simulation_id   '_fig003'],'-depsc2','-painters')

    %% 4.4) Plot posterior densities

    figure(004); clf
    mcmcplot(chain,[],res,'denspanel');
    saveas(gcf, [filepath_save_figs simulation_id   '_fig004'  '.fig'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig004'  '.pdf'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig004'  '.png'])
    print([filepath_save_figs simulation_id   '_fig004'],'-depsc2','-painters')

    figure(005); clf;
    mcmcplot_density(chain,res,filepath_save_figs,simulation_id);
    saveas(gcf, [filepath_save_figs simulation_id  '_fig005' '.fig'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig005' '.pdf'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig005' '.png'])
    print([filepath_save_figs simulation_id '_fig005'],'-depsc2','-painters')


    %% 4.5) Plot the MCMC pairs

    figure(006); clf
    mcmcplot(chain,[],res,'pairs');
    saveas(gcf, [filepath_save_figs simulation_id   '_fig006'  '.fig'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig006'  '.pdf'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig006'  '.png'])
    print([filepath_save_figs  '_fig006'],'-depsc2','-painters')

    %% 4.6) Plot the mcmc chain error standard deviation

    figure(007); clf
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
    saveas(gcf, [filepath_save_figs   simulation_id '_fig007'  '.fig'])
    saveas(gcf, [filepath_save_figs   simulation_id '_fig007'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig007'  '.png'])
    print([filepath_save_figs simulation_id '_fig007'],'-depsc2','-painters')

    %% 4.7) Plot the experimental data v the model simulated with the mean of the chain

    x = linspace(2,10)';
    figure(008)
    plot(data.xdata,data.ydata,'s');
    hold on
    plot(x,modelfun(x,mean(chain)),'-k')
    hold off
    legend({'data','model'},'Location','best')
    saveas(gcf, [filepath_save_figs simulation_id   '_fig008'  '.fig'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig008'  '.pdf'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig008'  '.png'])
    print([filepath_save_figs  '_fig008'],'-depsc2','-painters')

    %%  4.8) Plot the prediction intervals vs experimental data (scatter plot)

    figure(009); clf
    out = mcmcpred(res,chain,s2chain,x,modelfun);
    mcmcpredplot(out);
    hold on
    plot(data.xdata,data.ydata,'s'); % add data points to the plot
    xlim([0,10])
    xlabel('time [days]')
    ylabel('Rd [micrometres]')
    hold off
    title('Predictive envelopes of the model')
    ylim([0,500])
    box on
    grid on
    xticks([0,2,4,6,8,10])
    yticks([0,100,200,300,400,500])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig009'  '.fig'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig009'  '.pdf'])
    saveas(gcf, [filepath_save_figs simulation_id   '_fig009'  '.png'])
    print([filepath_save_figs simulation_id  '_fig009'],'-depsc2','-painters')


end

% change directory back to 2Code
cd ../
cd 2Code