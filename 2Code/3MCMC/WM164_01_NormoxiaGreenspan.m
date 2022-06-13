%% 2022--Murphy et al--Growth and adaptation mechanisms of tumour spheroids with time-dependent oxygen availability
% NormoxiaGreenspan

% This script uses the GitHub package developed by Marko Laine (https://mjlaine.github.io/mcmcstat/)

%% CODE STRUCTURE

% 0) Filepath to save results
% 1) Load data
% 2) Estimate initial outer radius, Ro(0)
% 3) Estimate outer radius when necrotic region forms, Rc
% 4) Estimate outer radius when inhibited region forms, \mathcal{R},
% referred to as Rd
% 6) Estimate all parameters (Ro(0), Rc,\mathcal{R}, s, lambda)
% 6.1) Global minimisation for first guess
% 6.2) Compute mse of model simulated with global min
% 6.3) Run MCMC chains
% 6.4) MCMC Rhat, mean and quantiles of chain
% 6.5) Sample MCMC chain to plot prediction intervals
% 6.6) Estimate mean and standard deviation to MCMC chain parameters to inform parameter bounds in deoxygenation, re-oxygenation simulations
% 7) Save data to .mat file
% 8) Plots
% 8.0) Load data set from 7
% 8.1) Plot fit to initial outer radius
% 8.2) Plot fit to Rc
% 8.3) Plot fit to Rd (\mathcalR)
% 8.4) plot the experimental data
% 8.5) plot the model simulated with the global min vs the experimental data
% 8.6) Plot chains
% 8.7) Plot posterior densities
% 8.8) Plot MCMC pairs
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
if ~exist('3MATLABsavefiles\3MCMC\WM164\1Normoxia', 'dir')
    mkdir(['3MATLABsavefiles\3MCMC\WM164\1Normoxia']);
end
if ~exist('4Figures\3MCMC\WM164\1Normoxia', 'dir')
    mkdir(['4Figures\3MCMC\WM164\1Normoxia']);
end


filepath_save = [newdir '\3MATLABsavefiles\3MCMC\WM164\1Normoxia\'];
filepath_save_figs =  [newdir '\4Figures\3MCMC\WM164\1Normoxia\'];

simulation_id = 'WM164_01_NormoxiaGreenspan';

% open directory with the MCMC functions
cd 'mcmcstat-master'

%% 1) Load data

% open data from 1Data folder
filepath_load = [newdir '\3MATLABsavefiles\1FirstAnalysisExperimentalData\WM164\'];

% remove outliers
WM164_outer_tmp = load([filepath_load 'WM164_confocal_scenario_1_outer.mat']);
WM164_necrotic_tmp = load([filepath_load 'WM164_confocal_scenario_1_necrotic.mat']);
WM164_inhibited_tmp = load([filepath_load 'WM164_confocal_scenario_1_inhibited.mat']);

WM164_outer.Day = WM164_outer_tmp.Day(WM164_outer_tmp.Outlier == 0);
WM164_outer.Radius = WM164_outer_tmp.Radius(WM164_outer_tmp.Outlier == 0);
WM164_outer.Outlier = WM164_outer_tmp.Outlier(WM164_outer_tmp.Outlier == 0);

WM164_necrotic.Day = WM164_necrotic_tmp.Day(WM164_necrotic_tmp.Outlier == 0);
WM164_necrotic.Radius = WM164_necrotic_tmp.Radius(WM164_necrotic_tmp.Outlier == 0);
WM164_necrotic.Outlier = WM164_necrotic_tmp.Outlier(WM164_necrotic_tmp.Outlier == 0);

WM164_inhibited.Day = WM164_inhibited_tmp.Day(WM164_inhibited_tmp.Outlier == 0);
WM164_inhibited.Radius = WM164_inhibited_tmp.Radius(WM164_inhibited_tmp.Outlier == 0);
WM164_inhibited.Outlier = WM164_inhibited_tmp.Outlier(WM164_inhibited_tmp.Outlier == 0);

%% 2) Estimate initial outer radius, Ro(0)

% obtain first day
day_min = min(WM164_outer.Day);

% select all outer radii data on the first day (and exclude outliers)
OuterRadius_initial = WM164_outer.Radius(logical((WM164_outer.Outlier==0).*(WM164_outer.Day==day_min)));

% fit a normal distribution to data
OuterRadiusInitial_normaldistribution_fit = fitdist(OuterRadius_initial,'Normal');

% output mean and standard deviation
OuterRadiusInitial_normaldistribution_fit.mu;
OuterRadiusInitial_normaldistribution_fit.sigma;

%% 3) Estimate outer radius when necrotic region forms, Rc

num_spheroids = length(WM164_outer.Radius);

Rcrit_estimate_per_spheroid = [];
Rcrit_estimate_counter = 1;
for i=1:num_spheroids
    if (WM164_outer_tmp.Outlier(i)==0) && (WM164_necrotic_tmp.Outlier(i)==0) % include if measurement is not an outlier
        OuterRadius_this_loop = WM164_outer_tmp.Radius(i);
        NecroticRadius_this_loop = WM164_necrotic_tmp.Radius(i);

        % for each spheroid compute R_c [micrometres]
        Rcrit_estimate_per_spheroid(Rcrit_estimate_counter) = sqrt((OuterRadius_this_loop.^2 - NecroticRadius_this_loop.^2) - (2*NecroticRadius_this_loop.^2/OuterRadius_this_loop)*(OuterRadius_this_loop - NecroticRadius_this_loop));

        Rcrit_estimate_counter = Rcrit_estimate_counter + 1;
    end
end

Rcrit_estimate_per_spheroid=Rcrit_estimate_per_spheroid';

% fit a normal distribution to data
Rc_normaldistribution_fit = fitdist(Rcrit_estimate_per_spheroid,'Normal');

% output mean and standard deviation
Rc_normaldistribution_fit.mu;
Rc_normaldistribution_fit.sigma;

%% 4) Estimate outer radius when inhibited region forms, \mathcal{R} (refer to as Rd)

num_spheroids = length(WM164_outer.Radius);

Rd_estimate_per_spheroid = zeros(1,1);
Rd_counter = 1;
for i=1:num_spheroids

    if (WM164_outer_tmp.Outlier(i)==0) && (WM164_necrotic_tmp.Outlier(i)==0) && (WM164_inhibited_tmp.Outlier(i)==0) % include if measurement is not an outlier

        OuterRadius_this_loop = WM164_outer_tmp.Radius(i);
        NecroticRadius_this_loop = WM164_necrotic_tmp.Radius(i);
        InhibitedRadius_this_loop = WM164_inhibited_tmp.Radius(i);

        %for each spheroid compute R_d [micrometres]
        if InhibitedRadius_this_loop > 0 % if in phase 2 or 3
            Rd_estimate_per_spheroid(Rd_counter) = sqrt((OuterRadius_this_loop.^2 - InhibitedRadius_this_loop.^2) - 2.*NecroticRadius_this_loop.^3./InhibitedRadius_this_loop + 2.*NecroticRadius_this_loop.^3./OuterRadius_this_loop);
            Rd_counter = Rd_counter + 1;
        end
    end
end
Rd_estimate_per_spheroid =Rd_estimate_per_spheroid';

% fit a normal distribution to data
Rd_normaldistribution_fit = fitdist(Rd_estimate_per_spheroid,'Normal');

% output mean and standard deviation
Rd_normaldistribution_fit.mu;
Rd_normaldistribution_fit.sigma;


%% 6) Estimate all parameters with MCMC (Ro(0), Rc,\mathcal{R}, s, lambda)

% theta(1) = Rd;
% theta(2) = lambda;
% theta(3) = s;
% theta(4) = Rc;
% theta(5) = R_o_init;

day_shift = 2;

data.xdata = [WM164_outer.Day;WM164_necrotic.Day;WM164_inhibited.Day]-day_shift; % Days
data.ydata = [WM164_outer.Radius;WM164_necrotic.Radius;WM164_inhibited.Radius]; % Outer, Necrotic, Inhibited
variable_id_vec = [ones(length(WM164_outer.Day),1);2*ones(length(WM164_necrotic.Day),1);3*ones(length(WM164_inhibited.Day),1) ];


modelfun = @(x,theta) function_Greenspan_1simulation(theta,...
    x,...
    variable_id_vec,...
    max(x));

ssfun    = @(theta,data) sum((data.ydata-modelfun(data.xdata,theta)).^2);

%% 6.1) Global minimisation for first guess

theta_first_guess =  [Rd_normaldistribution_fit.mu; % Rd
    0.1; % lambda
    0.15; % s
    Rc_normaldistribution_fit.mu;  % Rc
    OuterRadiusInitial_normaldistribution_fit.mu]; % Roinit

fun_likelihood = @(p) ssfun(p,data);
p_first_guess = theta_first_guess;
p_lower_bounds = [max(0,Rd_normaldistribution_fit.mu - 5*Rd_normaldistribution_fit.sigma); % theta(1) = Rd
    0.001; % theta(2) = lambda
    0.001; % theta(3) = s
    max(0,Rc_normaldistribution_fit.mu - 5*Rc_normaldistribution_fit.sigma);  % theta(4) = Rc
    max(0,OuterRadiusInitial_normaldistribution_fit.mu - 5*OuterRadiusInitial_normaldistribution_fit.sigma)]; % theta(5) = Roinit
p_upper_bounds =  [max(0,Rd_normaldistribution_fit.mu + 5*Rd_normaldistribution_fit.sigma); % theta(1) = Rd
    6; % theta(2) = lambda
    2; % theta(3) = s
    max(0,Rc_normaldistribution_fit.mu +5*Rc_normaldistribution_fit.sigma);  % theta(4) = Rc
    max(0,OuterRadiusInitial_normaldistribution_fit.mu + 5*OuterRadiusInitial_normaldistribution_fit.sigma)]; % theta(5) = Roinit

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

x_data_outer_pos = 1:length(WM164_outer.Day);
x_data_necrotic_pos = (length(WM164_outer.Day)+1):(length(WM164_outer.Day)+length(WM164_necrotic.Day));
x_data_inhibited_pos = (length(WM164_outer.Day)+length(WM164_necrotic.Day)+1):(length(WM164_outer.Day)+length(WM164_necrotic.Day)+length(WM164_inhibited.Day));

xdata_mse = [data.xdata(x_data_outer_pos);data.xdata(x_data_necrotic_pos);data.xdata(x_data_inhibited_pos)];
variable_id_vec_mse = [ones(length(data.xdata(x_data_outer_pos)),1);2*ones(length(data.xdata(x_data_necrotic_pos)),1);3*ones(length(data.xdata(x_data_inhibited_pos)),1) ];

[ydata_mse] = function_Greenspan_1simulation(p_mle,...
    [xdata_mse,xdata_mse,xdata_mse],...
    variable_id_vec_mse,...
    max(xdata_mse));

count_for_mse = length(data.xdata);

% MSE = (1/n)(sumi=1..N (yi-Yi)^2) where yi are observations and Yi is prediction
mse = (1./count_for_mse).*sum((data.ydata - ydata_mse).^2);
save([filepath_save 'MCMC_' simulation_id '_globalmin' '.mat'],'-v7.3');

%% 6.3) Run MCMC

params = {
    {'Rd',p_mle(1), p_lower_bounds(1), p_upper_bounds(1),Rd_normaldistribution_fit.mu,3*Rd_normaldistribution_fit.sigma}
    {'lambda', p_mle(2), p_lower_bounds(2), p_upper_bounds(2)}
    {'s', p_mle(3), p_lower_bounds(3), p_upper_bounds(3)}
    {'Rc',p_mle(4),p_lower_bounds(4), p_upper_bounds(4),Rc_normaldistribution_fit.mu,3*Rc_normaldistribution_fit.sigma}
    {'Roinit',p_mle(5),p_lower_bounds(5), p_upper_bounds(5),OuterRadiusInitial_normaldistribution_fit.mu,3*OuterRadiusInitial_normaldistribution_fit.sigma}
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

time_prediction_interval_one_variable = linspace(0,10)';
len_time_prediction_interval_one_variable = length(time_prediction_interval_one_variable);
time_prediction_interval = [time_prediction_interval_one_variable;time_prediction_interval_one_variable;time_prediction_interval_one_variable];

x_data_outer_pos_prediction_interval = 1:length(time_prediction_interval_one_variable);
x_data_necrotic_pos_prediction_interval = (length(time_prediction_interval_one_variable)+1):(length(time_prediction_interval_one_variable)+length(time_prediction_interval_one_variable));
x_data_inhibited_pos_prediction_interval = (length(time_prediction_interval_one_variable)+length(time_prediction_interval_one_variable)+1):(length(time_prediction_interval_one_variable)+length(time_prediction_interval_one_variable)+length(time_prediction_interval_one_variable));

variable_id_vec_prediction_interval = [ones(len_time_prediction_interval_one_variable,1);2*ones(len_time_prediction_interval_one_variable,1);3*ones(len_time_prediction_interval_one_variable,1) ];

modelfun_prediction_interval = @(x,theta) function_Greenspan_1simulation(theta,...
    x,...
    variable_id_vec_prediction_interval,...
    max(x));

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

%% 6.6) Esimate mean and standard deviation to MCMC chain parameters to inform parameter bounds in deoxygenation, re-oxygenation simulations

params_withnormalfit = params;

for r=1:size(params_withnormalfit,1)
    params_withnormalfit{r}{7} = chain_mean(r);
    params_withnormalfit{r}{8} = chain_std(r);
end


%% 7) Save outputs

close all
save([filepath_save 'MCMC_WM164_01_NormoxiaGreenspan' '.mat'],'-v7.3');


%% 8) PLOTS

if include_plots == 1

    %% 8.0) Load experimental data set

    load([filepath_save 'MCMC_WM164_01_NormoxiaGreenspan' '.mat']);


    %% 8.1) Plot fit to initial outer radius

    figure(001); clf;
    hold on
    OuterRadius_initial_plot = linspace(OuterRadiusInitial_normaldistribution_fit.mu-3.*OuterRadiusInitial_normaldistribution_fit.sigma,...
        OuterRadiusInitial_normaldistribution_fit.mu+3.*OuterRadiusInitial_normaldistribution_fit.sigma);
    y_OuterRadiusInitial_normaldistribution_fit = pdf(OuterRadiusInitial_normaldistribution_fit,OuterRadius_initial_plot);
    plot(OuterRadius_initial_plot,y_OuterRadiusInitial_normaldistribution_fit,'r')
    xlim([0,500])
    xlabel('Roinit')
    ylabel(['probability'])
    xticks([0,100,200,300,400,500])
    box on
    grid on
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig001'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig001'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig001'  '.png'])
    print([filepath_save_figs simulation_id '_fig001'],'-depsc2','-painters')

    %% 8.2) Plot fit to Rc

    figure(002); clf;
    hold on
    Rc_plot = linspace(Rc_normaldistribution_fit.mu-3.*Rc_normaldistribution_fit.sigma,...
        Rc_normaldistribution_fit.mu+3.*Rc_normaldistribution_fit.sigma);
    y_Rc_normaldistribution_fit = pdf(Rc_normaldistribution_fit,Rc_plot);
    plot(Rc_plot,y_Rc_normaldistribution_fit,'r')
    xlim([0,500])
    xlabel('Rc')
    ylabel(['probability'])
    box on
    xticks([0,100,200,300,400,500])
    grid on
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig002'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig002'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig002'  '.png'])
    print([filepath_save_figs simulation_id '_fig002'],'-depsc2','-painters')

    %% 8.3) Plot fit to Rd

    figure(003); clf;
    hold on
    Rd_plot = linspace(Rd_normaldistribution_fit.mu-3.*Rd_normaldistribution_fit.sigma,...
        Rd_normaldistribution_fit.mu+3.*Rd_normaldistribution_fit.sigma);
    y_Rd_normaldistribution_fit = pdf(Rd_normaldistribution_fit,Rd_plot);
    plot(Rd_plot,y_Rd_normaldistribution_fit,'r')
    xlim([0,500])
    xlabel('Rd')
    ylabel(['probability'])
    box on
    xticks([0,100,200,300,400,500])
    grid on
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig003'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig003'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig003'  '.png'])
    print([filepath_save_figs simulation_id '_fig003'],'-depsc2','-painters')



    %% 8.4) Plot experimental data

    figure(004); clf
    hold on
    scatter(data.xdata(x_data_outer_pos)+day_shift,data.ydata(x_data_outer_pos),'g')
    scatter(data.xdata(x_data_necrotic_pos)+day_shift,data.ydata(x_data_necrotic_pos),'k')
    scatter(data.xdata(x_data_inhibited_pos)+day_shift,data.ydata(x_data_inhibited_pos),'m')
    xlim([0,10])
    ylim([0,500])
    xlabel('time [days]')
    ylabel(['Radius'])
    box on
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig004'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig004'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig004'  '.png'])
    print([filepath_save_figs simulation_id '_fig004'],'-depsc2','-painters')

    %% 8.5) plot the model simulated with the global min vs the experimental data

    xdata_globalmin = 0:0.1:10;
    variable_id_vec_globalmin = [ones(length(xdata_globalmin),1);2*ones(length(xdata_globalmin),1);3*ones(length(xdata_globalmin),1) ];
    x_data_outer_pos_globalmin = 1:length(xdata_globalmin);
    x_data_necrotic_pos_globalmin = (length(xdata_globalmin)+1):(length(xdata_globalmin)+length(xdata_globalmin));
    x_data_inhibited_pos_globalmin = (length(xdata_globalmin)+length(xdata_globalmin)+1):(length(xdata_globalmin)+length(xdata_globalmin)+length(xdata_globalmin));

    [ydata_globalmin] = function_Greenspan_1simulation(p_mle,...
        [xdata_globalmin,xdata_globalmin,xdata_globalmin],...
        variable_id_vec_globalmin,...
        max(xdata_globalmin));

    figure(005); clf
    hold on
    scatter(data.xdata(x_data_outer_pos)+day_shift,data.ydata(x_data_outer_pos),'g')
    scatter(data.xdata(x_data_necrotic_pos)+day_shift,data.ydata(x_data_necrotic_pos),'k')
    scatter(data.xdata(x_data_inhibited_pos)+day_shift,data.ydata(x_data_inhibited_pos),'m')
    hold on
    plot(xdata_globalmin+day_shift,ydata_globalmin(x_data_outer_pos_globalmin),'g--')
    plot(xdata_globalmin+day_shift,ydata_globalmin(x_data_necrotic_pos_globalmin),'k--')
    plot(xdata_globalmin+day_shift,ydata_globalmin(x_data_inhibited_pos_globalmin),'m--')
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

    %% 8.8) Plot MCMC pairs

    figure(009); clf
    mcmcplot(chain,[],res,'pairs');
    saveas(gcf, [filepath_save_figs   simulation_id '_fig009'  '.fig'])
    saveas(gcf, [filepath_save_figs   simulation_id '_fig009'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig009'  '.png'])
    print([filepath_save_figs simulation_id '_fig009'],'-depsc2','-painters')

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

    xdata_meanchain = 0:0.1:10;
    variable_id_vec_meanchain = [ones(length(xdata_meanchain),1);2*ones(length(xdata_meanchain),1);3*ones(length(xdata_meanchain),1) ];

    [ydata_meanchain] = function_Greenspan_1simulation(mean(chain),...
        [xdata_meanchain,xdata_meanchain,xdata_meanchain],...
        variable_id_vec_meanchain,...
        max(xdata_meanchain));

    x_data_outer_pos_meanchain = 1:length(xdata_meanchain);
    x_data_necrotic_pos_meanchain = (length(xdata_meanchain)+1):(length(xdata_meanchain)+length(xdata_meanchain));
    x_data_inhibited_pos_meanchain = (length(xdata_meanchain)+length(xdata_meanchain)+1):(length(xdata_meanchain)+length(xdata_meanchain)+length(xdata_meanchain));

    figure(011); clf
    hold on
    scatter(data.xdata(x_data_outer_pos)+day_shift,data.ydata(x_data_outer_pos),'g')
    scatter(data.xdata(x_data_necrotic_pos)+day_shift,data.ydata(x_data_necrotic_pos),'k')
    scatter(data.xdata(x_data_inhibited_pos)+day_shift,data.ydata(x_data_inhibited_pos),'m')
    hold on
    plot(xdata_meanchain+day_shift,ydata_meanchain(x_data_outer_pos_meanchain),'g--')
    plot(xdata_meanchain+day_shift,ydata_meanchain(x_data_necrotic_pos_meanchain),'k--')
    plot(xdata_meanchain+day_shift,ydata_meanchain(x_data_inhibited_pos_meanchain),'m--')
    xlim([0,10])
    ylim([0,500])
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
    fillyy(time_prediction_interval(x_data_outer_pos_prediction_interval)+day_shift,predlims(1,x_data_outer_pos_prediction_interval),predlims(2*nn-1,x_data_outer_pos_prediction_interval),dimc_outer);
    hold on
    for k=2:(nn-1)
        fillyy(time_prediction_interval(x_data_outer_pos_prediction_interval)+day_shift,predlims(k,x_data_outer_pos_prediction_interval),predlims(2*nn-k,x_data_outer_pos_prediction_interval),dimc_outer.*0.9.^(k-1));
    end

    % INHIBITED
    dimc_inhibited = [1.0 0.8 0.9]; % dimmest (lightest) color INHIBITED - magenta
    fillyy(time_prediction_interval(x_data_inhibited_pos_prediction_interval)+day_shift,predlims(1,x_data_inhibited_pos_prediction_interval),predlims(2*nn-1,x_data_inhibited_pos_prediction_interval),dimc_inhibited);
    hold on
    for k=2:(nn-1)
        fillyy(time_prediction_interval(x_data_inhibited_pos_prediction_interval)+day_shift,predlims(k,x_data_inhibited_pos_prediction_interval),predlims(2*nn-k,x_data_inhibited_pos_prediction_interval),dimc_inhibited.*0.9.^(k-1));
    end

    % NECROTIC
    dimc_necrotic = [0.8 0.8 0.8]; % dimmest (lightest) color NECROTIC - grey
    fillyy(time_prediction_interval(x_data_necrotic_pos_prediction_interval)+day_shift,predlims(1,x_data_necrotic_pos_prediction_interval),predlims(2*nn-1,x_data_necrotic_pos_prediction_interval),dimc_necrotic);
    hold on
    for k=2:(nn-1)
        fillyy(time_prediction_interval(x_data_necrotic_pos_prediction_interval)+day_shift,predlims(k,x_data_necrotic_pos_prediction_interval),predlims(2*nn-k,x_data_necrotic_pos_prediction_interval),dimc_necrotic.*0.9.^(k-1));
    end

    % PLOT experimental data
    hold on
    scatter(data.xdata(x_data_outer_pos)+day_shift,data.ydata(x_data_outer_pos),'g')
    scatter(data.xdata(x_data_necrotic_pos)+day_shift,data.ydata(x_data_necrotic_pos),'k')
    scatter(data.xdata(x_data_inhibited_pos)+day_shift,data.ydata(x_data_inhibited_pos),'m')

    % figure properties
    xlim([0,10])
    ylim([0,500])
    box on
    xlabel('time [days]')
    ylabel('Radius [microns]')
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig012'  '.fig'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig012'  '.pdf'])
    saveas(gcf, [filepath_save_figs  simulation_id  '_fig012'  '.png'])
    print([filepath_save_figs simulation_id '_fig012'],'-depsc2','-painters')

    %%  8.12) Plot the prediction intervals vs experimental data (box charts)

    fig_boxchart_predinterval = openfig([newdir '\4Figures\1FirstAnalysisExperimentalData\WM164\Scenario_1_OIN_1.fig']);

    figure(fig_boxchart_predinterval)
    hold on
    np = size(predlims,1);
    nn = (np+1)/2; % median
    np = nn-1;
    % OUTER
    dimc_outer = [0.8 1.0 0.8]; % dimmest (lightest) color OUTER- green
    fillyy(time_prediction_interval(x_data_outer_pos_prediction_interval)+day_shift,predlims(1,x_data_outer_pos_prediction_interval),predlims(2*nn-1,x_data_outer_pos_prediction_interval),dimc_outer);
    hold on
    for k=2:(nn-1)
        fillyy(time_prediction_interval(x_data_outer_pos_prediction_interval)+day_shift,predlims(k,x_data_outer_pos_prediction_interval),predlims(2*nn-k,x_data_outer_pos_prediction_interval),dimc_outer.*0.9.^(k-1));
    end

    % INHIBITED
    dimc_inhibited = [1.0 0.8 0.9]; % dimmest (lightest) color INHIBITED - magenta
    fillyy(time_prediction_interval(x_data_inhibited_pos_prediction_interval)+day_shift,predlims(1,x_data_inhibited_pos_prediction_interval),predlims(2*nn-1,x_data_inhibited_pos_prediction_interval),dimc_inhibited);
    hold on
    for k=2:(nn-1)
        fillyy(time_prediction_interval(x_data_inhibited_pos_prediction_interval)+day_shift,predlims(k,x_data_inhibited_pos_prediction_interval),predlims(2*nn-k,x_data_inhibited_pos_prediction_interval),dimc_inhibited.*0.9.^(k-1));
    end

    % NECROTIC
    dimc_necrotic = [0.8 0.8 0.8]; % dimmest (lightest) color NECROTIC - grey
    fillyy(time_prediction_interval(x_data_necrotic_pos_prediction_interval)+day_shift,predlims(1,x_data_necrotic_pos_prediction_interval),predlims(2*nn-1,x_data_necrotic_pos_prediction_interval),dimc_necrotic);
    hold on
    for k=2:(nn-1)
        fillyy(time_prediction_interval(x_data_necrotic_pos_prediction_interval)+day_shift,predlims(k,x_data_necrotic_pos_prediction_interval),predlims(2*nn-k,x_data_necrotic_pos_prediction_interval),dimc_necrotic.*0.9.^(k-1));
    end


    % figure properties
    xlim([0,10])
    ylim([0,500])
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
