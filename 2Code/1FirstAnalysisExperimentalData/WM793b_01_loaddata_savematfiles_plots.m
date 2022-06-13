%% Convert confocal WM793b data to .mat file and plot experimental data
% 0.1) Load data
% 1) PLOTS - Scenario 1 - 21% all days
% 2) PLOTS - Scenario 2 - 2% all days
% 3) PLOTS - Scenario 3 - Deoxygenation - 21% to Day 2 then 2%
% 4) PLOTS - Scenario 4 - Re-oxygenation - 2% to Day 2 then 21%
% 5) PLOTS - Scenario 5 - Re-oxygenation - 2% to Day 4 then 21%
% 6) PLOTS - Compare Scenario 1 and 2
% 7) Save data to .mat files
% 8) Change directory back to 2Code

%% 0.1) Load WM793b data

% open data from 1Data folder
mydir  = pwd; idcs   = strfind(mydir,'\'); newdir = mydir(1:idcs(end)-1);

confocal.AllData = readtable([newdir '\1Data\ExpOA_WM793b.xlsx'],'Sheet', 'Sheet1'); % load data
wm793b_data_include = (sum(cell2mat(confocal.AllData.CellLine) == 'WM793',2)==5); % identify WM793b data
confocal_WM793b.AllDatatmp = confocal.AllData(wm793b_data_include,:); % load WM793b data
confocal_WM793b.AllData = confocal_WM793b.AllDatatmp((confocal_WM793b.AllDatatmp.Keep==1),:); % only keep those that are not multiple spheroids or damaged spheroids

% create folders
cd ../
if ~exist('3MATLABsavefiles\1FirstAnalysisExperimentalData\WM793b', 'dir')
    mkdir(['3MATLABsavefiles\1FirstAnalysisExperimentalData\WM793b']);
end
if ~exist('4Figures\1FirstAnalysisExperimentalData\WM793b', 'dir')
    mkdir(['4Figures\1FirstAnalysisExperimentalData\WM793b']);
end

cd '2Code';
cd '1FirstAnalysisExperimentalData';

filepath_save = [newdir '\3MATLABsavefiles\1FirstAnalysisExperimentalData\WM793b\'];
filepath_save_figs =  [newdir '\4Figures\1FirstAnalysisExperimentalData\WM793b\'];

% setting for scatter plots
jitter_plot = 0.18;


%% 1) PLOTS -  Scenario 1 - 21% all days

scenario_1_conditions = [17,18,19,15,16];


[scenario_1_plotdata.OuterRadius_to_plot,...
    scenario_1_plotdata.InhibitedRadius_to_plot,...
    scenario_1_plotdata.NecroticRadius_to_plot,...
    scenario_1_plotdata.PIMOuterRadius_to_plot,...
    scenario_1_plotdata.PIMInnerRadius_to_plot,...
    scenario_1_plotdata.Days_to_plot,...
    scenario_1_plotdata.OuterRadius_mean_to_plot,...
    scenario_1_plotdata.InhibitedRadius_mean_to_plot,...
    scenario_1_plotdata.NecroticRadius_mean_to_plot,...
    scenario_1_plotdata.PIMOuterRadius_mean_to_plot,...
    scenario_1_plotdata.PIMInnerRadius_mean_to_plot,...
    scenario_1_plotdata.Days_unique_to_plot,...
    scenario_1_plotdata.OuterRadius_std_to_plot,...
    scenario_1_plotdata.InhibitedRadius_std_to_plot,...
    scenario_1_plotdata.NecroticRadius_std_to_plot,...
    scenario_1_plotdata.PIMOuterRadius_std_to_plot,...
    scenario_1_plotdata.PIMInnerRadius_std_to_plot,...
    scenario_1_plotdata.OuterRadius_count_to_plot,...
    scenario_1_plotdata.InhibitedRadius_count_to_plot,...
    scenario_1_plotdata.NecroticRadius_count_to_plot,...
    scenario_1_plotdata.PIMOuterRadius_count_to_plot,...
    scenario_1_plotdata.PIMInnerRadius_count_to_plot]  =  function_plotdata_scenario(scenario_1_conditions,...
    confocal_WM793b.AllData.Condition,...
    confocal_WM793b.AllData.OuterRadius,...
    confocal_WM793b.AllData.InhibitedRadius,...
    confocal_WM793b.AllData.NecroticRadius,...
    confocal_WM793b.AllData.PIMOuterRadius,...
    confocal_WM793b.AllData.PIMInnerRadius,...
    confocal_WM793b.AllData.Day);


% BOX PLOTS
figure
hold on
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 2),scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 3),scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 4),scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 6),scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 8),scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
hold on
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 2),scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 3),scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 4),scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 6),scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 8),scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
hold on
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 2),scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 3),scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 4),scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 6),scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 8),scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')

% exclude outliers from scatter
%use isoutlier function
scenario_1_days_to_plot_outer_tmp =scenario_1_plotdata.Days_to_plot;
scenario_1_plotdata_outer_tmp =  scenario_1_plotdata.OuterRadius_to_plot;
scenario_1_days_to_plot_inhibited_tmp =scenario_1_plotdata.Days_to_plot;
scenario_1_plotdata_inhibited_tmp = scenario_1_plotdata.InhibitedRadius_to_plot;
scenario_1_days_to_plot_necrotic_tmp = scenario_1_plotdata.Days_to_plot;
scenario_1_plotdata_necrotic_tmp = scenario_1_plotdata.NecroticRadius_to_plot;

[scenario_1_days_to_plot_outer,scenario_1_plotdata_outer] = function_removeoutliers_for_plot(scenario_1_days_to_plot_outer_tmp,scenario_1_plotdata_outer_tmp);
[scenario_1_days_to_plot_inhibited,scenario_1_plotdata_inhibited] = function_removeoutliers_for_plot(scenario_1_days_to_plot_inhibited_tmp,scenario_1_plotdata_inhibited_tmp);
[scenario_1_days_to_plot_necrotic,scenario_1_plotdata_necrotic] = function_removeoutliers_for_plot(scenario_1_days_to_plot_necrotic_tmp,scenario_1_plotdata_necrotic_tmp);

s3 = scatter(scenario_1_days_to_plot_outer,scenario_1_plotdata_outer,'g','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s2 = scatter(scenario_1_days_to_plot_inhibited,scenario_1_plotdata_inhibited,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(scenario_1_days_to_plot_necrotic,scenario_1_plotdata_necrotic,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);

alpha(s3,0.2)
alpha(s2,0.2)
alpha(s1,0.2)
xlim([0,8.5])
ylim([0,400])
xlabel('time [days]')
ylabel('Radius [microns]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_1_OIN_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_1_OIN_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_1_OIN_1'  '.png'])
print([filepath_save_figs  'Scenario_1_OIN_1'],'-depsc2','-painters')

% BOX PLOT WITH PIM

figure
hold on
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 2),scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 3),scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 4),scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 6),scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 8),scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
hold on
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 2),scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 3),scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 4),scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 6),scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 8),scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
hold on
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 2),scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 3),scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 4),scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 6),scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 8),scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
hold on
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 2),scenario_1_plotdata.PIMOuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 3),scenario_1_plotdata.PIMOuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 4),scenario_1_plotdata.PIMOuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 6),scenario_1_plotdata.PIMOuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 8),scenario_1_plotdata.PIMOuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])

% exclude outliers from scatter
%use isoutlier function
scenario_1_days_to_plot_outer_tmp =scenario_1_plotdata.Days_to_plot;
scenario_1_plotdata_outer_tmp =  scenario_1_plotdata.OuterRadius_to_plot;
scenario_1_days_to_plot_inhibited_tmp =scenario_1_plotdata.Days_to_plot;
scenario_1_plotdata_inhibited_tmp = scenario_1_plotdata.InhibitedRadius_to_plot;
scenario_1_days_to_plot_necrotic_tmp = scenario_1_plotdata.Days_to_plot;
scenario_1_plotdata_necrotic_tmp = scenario_1_plotdata.NecroticRadius_to_plot;
scenario_1_days_to_plot_PIMOuter_tmp = scenario_1_plotdata.Days_to_plot;
scenario_1_plotdata_PIMOuter_tmp = scenario_1_plotdata.PIMOuterRadius_to_plot;

[scenario_1_days_to_plot_outer,scenario_1_plotdata_outer] = function_removeoutliers_for_plot(scenario_1_days_to_plot_outer_tmp,scenario_1_plotdata_outer_tmp);
[scenario_1_days_to_plot_inhibited,scenario_1_plotdata_inhibited] = function_removeoutliers_for_plot(scenario_1_days_to_plot_inhibited_tmp,scenario_1_plotdata_inhibited_tmp);
[scenario_1_days_to_plot_necrotic,scenario_1_plotdata_necrotic] = function_removeoutliers_for_plot(scenario_1_days_to_plot_necrotic_tmp,scenario_1_plotdata_necrotic_tmp);
[scenario_1_days_to_plot_PIMOuter,scenario_1_plotdata_PIMOuter] = function_removeoutliers_for_plot(scenario_1_days_to_plot_PIMOuter_tmp,scenario_1_plotdata_PIMOuter_tmp);

s3 = scatter(scenario_1_days_to_plot_outer,scenario_1_plotdata_outer,'g','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s2 = scatter(scenario_1_days_to_plot_inhibited,scenario_1_plotdata_inhibited,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(scenario_1_days_to_plot_necrotic,scenario_1_plotdata_necrotic,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s0 = scatter(scenario_1_days_to_plot_PIMOuter,scenario_1_plotdata_PIMOuter,'filled','MarkerFaceColor',[0,1,1],'jitter', 'on', 'jitterAmount', jitter_plot);

alpha(s3,0.2)
alpha(s2,0.2)
alpha(s1,0.2)
alpha(s0,0.2)
xlim([0,8.5])
ylim([0,400])
xlabel('time [days]')
ylabel('Radius [microns]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_1_OINP_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_1_OINP_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_1_OINP_1'  '.png'])
print([filepath_save_figs  'Scenario_1_OINP_1'],'-depsc2','-painters')

% BOX PLOTS for inhibited radius/outer radius, necrotic radius/outer radius,
% pim outer radius/inhibited radius

% BOX - inh/outer

figure
hold on
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 2),...
    scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 3),...
    scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 4),...
    scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 6),...
    scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 8),...
    scenario_1_plotdata.InhibitedRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')

        % exclude outliers from scatter
        %use isoutlier function
scenario_1_days_to_plot_inhibited_tmp =scenario_1_plotdata.Days_to_plot;
scenario_1_plotdata_inhdivout_tmp = scenario_1_plotdata.InhibitedRadius_to_plot./scenario_1_plotdata.OuterRadius_to_plot;
[scenario_1_days_to_plot_inhdivout,scenario_1_plotdata_inhdivout] = function_removeoutliers_for_plot(scenario_1_days_to_plot_inhibited_tmp,scenario_1_plotdata_inhdivout_tmp);
 s2 = scatter(scenario_1_days_to_plot_inhdivout,scenario_1_plotdata_inhdivout,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Inhibited Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_1_inhdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_1_inhdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_1_inhdivout_1'  '.png'])
print([filepath_save_figs  'Scenario_1_inhdivout_1'],'-depsc2','-painters')

% BOX - nec/outer

figure
hold on
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 2),...
    scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 3),...
    scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 4),...
    scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 6),...
    scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 8),...
    scenario_1_plotdata.NecroticRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')

        % exclude outliers from scatter
        %use isoutlier function
scenario_1_days_to_plot_necrotic_tmp =scenario_1_plotdata.Days_to_plot;
scenario_1_plotdata_necdivout_tmp = scenario_1_plotdata.NecroticRadius_to_plot./scenario_1_plotdata.OuterRadius_to_plot;
[scenario_1_days_to_plot_necdivout,scenario_1_plotdata_necdivout] = function_removeoutliers_for_plot(scenario_1_days_to_plot_necrotic_tmp,scenario_1_plotdata_necdivout_tmp);
 s2 = scatter(scenario_1_days_to_plot_necdivout,scenario_1_plotdata_necdivout,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Necrotic Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_1_necdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_1_necdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_1_necdivout_1'  '.png'])
print([filepath_save_figs  'Scenario_1_necdivout_1'],'-depsc2','-painters')


% BOX - PIM/outer

figure
hold on
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 2),...
    scenario_1_plotdata.PIMOuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 2),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 3),...
    scenario_1_plotdata.PIMOuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 3),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 4),...
    scenario_1_plotdata.PIMOuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 4),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 6),...
    scenario_1_plotdata.PIMOuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 6),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_1_plotdata.Days_to_plot(scenario_1_plotdata.Days_to_plot == 8),...
    scenario_1_plotdata.PIMOuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8)./scenario_1_plotdata.OuterRadius_to_plot(scenario_1_plotdata.Days_to_plot == 8),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])

        % exclude outliers from scatter
        %use isoutlier function
scenario_1_days_to_plot_pimouter_tmp =scenario_1_plotdata.Days_to_plot;
scenario_1_plotdata_pimouterdivout_tmp = scenario_1_plotdata.PIMOuterRadius_to_plot./scenario_1_plotdata.OuterRadius_to_plot;
[scenario_1_days_to_plot_pimouterdivout,scenario_1_plotdata_pimouterdivout] = function_removeoutliers_for_plot(scenario_1_days_to_plot_pimouter_tmp,scenario_1_plotdata_pimouterdivout_tmp);
 s2 = scatter(scenario_1_days_to_plot_pimouterdivout,scenario_1_plotdata_pimouterdivout,'filled','MarkerFaceColor',[0,1,1],'jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Hypoxic Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_1_pimdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_1_pimdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_1_pimdivout_1'  '.png'])
print([filepath_save_figs  'Scenario_1_pimdivout_1'],'-depsc2','-painters')


%% 2) PLOTS - Scenario 2 - 2% all days
scenario_2_conditions = [1,2,3,4];

[scenario_2_plotdata.OuterRadius_to_plot,...
    scenario_2_plotdata.InhibitedRadius_to_plot,...
    scenario_2_plotdata.NecroticRadius_to_plot,...
    scenario_2_plotdata.PIMOuterRadius_to_plot,...
    scenario_2_plotdata.PIMInnerRadius_to_plot,...
    scenario_2_plotdata.Days_to_plot,...
    scenario_2_plotdata.OuterRadius_mean_to_plot,...
    scenario_2_plotdata.InhibitedRadius_mean_to_plot,...
    scenario_2_plotdata.NecroticRadius_mean_to_plot,...
    scenario_2_plotdata.PIMOuterRadius_mean_to_plot,...
    scenario_2_plotdata.PIMInnerRadius_mean_to_plot,...
    scenario_2_plotdata.Days_unique_to_plot,...
    scenario_2_plotdata.OuterRadius_std_to_plot,...
    scenario_2_plotdata.InhibitedRadius_std_to_plot,...
    scenario_2_plotdata.NecroticRadius_std_to_plot,...
    scenario_2_plotdata.PIMOuterRadius_std_to_plot,...
    scenario_2_plotdata.PIMInnerRadius_std_to_plot,...
    scenario_2_plotdata.OuterRadius_count_to_plot,...
    scenario_2_plotdata.InhibitedRadius_count_to_plot,...
    scenario_2_plotdata.NecroticRadius_count_to_plot,...
    scenario_2_plotdata.PIMOuterRadius_count_to_plot,...
    scenario_2_plotdata.PIMInnerRadius_count_to_plot]  =  function_plotdata_scenario(scenario_2_conditions,...
    confocal_WM793b.AllData.Condition,...
    confocal_WM793b.AllData.OuterRadius,...
    confocal_WM793b.AllData.InhibitedRadius,...
    confocal_WM793b.AllData.NecroticRadius,...
    confocal_WM793b.AllData.PIMOuterRadius,...
    confocal_WM793b.AllData.PIMInnerRadius,...
    confocal_WM793b.AllData.Day);


% BOX PLOTS
figure
hold on
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 2),scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 3),scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 3),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 4),scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 6),scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 8),scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
hold on
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 2),scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 3),scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 4),scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 6),scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 8),scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
hold on
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 2),scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 3),scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 4),scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 6),scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 8),scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')

% exclude outliers from scatter
%use isoutlier function
scenario_2_days_to_plot_outer_tmp =scenario_2_plotdata.Days_to_plot;
scenario_2_plotdata_outer_tmp =  scenario_2_plotdata.OuterRadius_to_plot;
scenario_2_days_to_plot_inhibited_tmp =scenario_2_plotdata.Days_to_plot;
scenario_2_plotdata_inhibited_tmp = scenario_2_plotdata.InhibitedRadius_to_plot;
scenario_2_days_to_plot_necrotic_tmp = scenario_2_plotdata.Days_to_plot;
scenario_2_plotdata_necrotic_tmp = scenario_2_plotdata.NecroticRadius_to_plot;

[scenario_2_days_to_plot_outer,scenario_2_plotdata_outer] = function_removeoutliers_for_plot(scenario_2_days_to_plot_outer_tmp,scenario_2_plotdata_outer_tmp);
[scenario_2_days_to_plot_inhibited,scenario_2_plotdata_inhibited] = function_removeoutliers_for_plot(scenario_2_days_to_plot_inhibited_tmp,scenario_2_plotdata_inhibited_tmp);
[scenario_2_days_to_plot_necrotic,scenario_2_plotdata_necrotic] = function_removeoutliers_for_plot(scenario_2_days_to_plot_necrotic_tmp,scenario_2_plotdata_necrotic_tmp);

s3 = scatter(scenario_2_days_to_plot_outer,scenario_2_plotdata_outer,'g','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s2 = scatter(scenario_2_days_to_plot_inhibited,scenario_2_plotdata_inhibited,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(scenario_2_days_to_plot_necrotic,scenario_2_plotdata_necrotic,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);

alpha(s3,0.2)
alpha(s2,0.2)
alpha(s1,0.2)
xlim([0,8.5])
ylim([0,400])
xlabel('time [days]')
ylabel('Radius [microns]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_2_OIN_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_2_OIN_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_2_OIN_1'  '.png'])
print([filepath_save_figs  'Scenario_2_OIN_1'],'-depsc2','-painters')

% BOX PLOT WITH PIM

figure
hold on
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 2),scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 3),scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 3),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 4),scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 6),scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 8),scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
hold on
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 2),scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 3),scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 4),scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 6),scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 8),scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
hold on
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 2),scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 3),scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 4),scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 6),scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 8),scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
hold on
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 2),scenario_2_plotdata.PIMOuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 3),scenario_2_plotdata.PIMOuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 3),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 4),scenario_2_plotdata.PIMOuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 6),scenario_2_plotdata.PIMOuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 8),scenario_2_plotdata.PIMOuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])

% exclude outliers from scatter
%use isoutlier function
scenario_2_days_to_plot_outer_tmp =scenario_2_plotdata.Days_to_plot;
scenario_2_plotdata_outer_tmp =  scenario_2_plotdata.OuterRadius_to_plot;
scenario_2_days_to_plot_inhibited_tmp =scenario_2_plotdata.Days_to_plot;
scenario_2_plotdata_inhibited_tmp = scenario_2_plotdata.InhibitedRadius_to_plot;
scenario_2_days_to_plot_necrotic_tmp = scenario_2_plotdata.Days_to_plot;
scenario_2_plotdata_necrotic_tmp = scenario_2_plotdata.NecroticRadius_to_plot;
scenario_2_days_to_plot_PIMOuter_tmp = scenario_2_plotdata.Days_to_plot;
scenario_2_plotdata_PIMOuter_tmp = scenario_2_plotdata.PIMOuterRadius_to_plot;

[scenario_2_days_to_plot_outer,scenario_2_plotdata_outer] = function_removeoutliers_for_plot(scenario_2_days_to_plot_outer_tmp,scenario_2_plotdata_outer_tmp);
[scenario_2_days_to_plot_inhibited,scenario_2_plotdata_inhibited] = function_removeoutliers_for_plot(scenario_2_days_to_plot_inhibited_tmp,scenario_2_plotdata_inhibited_tmp);
[scenario_2_days_to_plot_necrotic,scenario_2_plotdata_necrotic] = function_removeoutliers_for_plot(scenario_2_days_to_plot_necrotic_tmp,scenario_2_plotdata_necrotic_tmp);
[scenario_2_days_to_plot_PIMOuter,scenario_2_plotdata_PIMOuter] = function_removeoutliers_for_plot(scenario_2_days_to_plot_PIMOuter_tmp,scenario_2_plotdata_PIMOuter_tmp);

s3 = scatter(scenario_2_days_to_plot_outer,scenario_2_plotdata_outer,'g','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s2 = scatter(scenario_2_days_to_plot_inhibited,scenario_2_plotdata_inhibited,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(scenario_2_days_to_plot_necrotic,scenario_2_plotdata_necrotic,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s0 = scatter(scenario_2_days_to_plot_PIMOuter,scenario_2_plotdata_PIMOuter,'filled','MarkerFaceColor',[0,1,1],'jitter', 'on', 'jitterAmount', jitter_plot);

alpha(s3,0.2)
alpha(s2,0.2)
alpha(s1,0.2)
alpha(s0,0.2)
xlim([0,8.5])
ylim([0,400])
xlabel('time [days]')
ylabel('Radius [microns]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_2_OINP_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_2_OINP_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_2_OINP_1'  '.png'])
print([filepath_save_figs  'Scenario_2_OINP_1'],'-depsc2','-painters')


% BOX PLOTS for inhibited radius/outer radius, necrotic radius/outer radius,
% pim outer radius/inhibited radius

% BOX - inh/outer

figure
hold on
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 2),...
    scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2)./scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 4),...
    scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4)./scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 6),...
    scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6)./scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 8),...
    scenario_2_plotdata.InhibitedRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8)./scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')

        % exclude outliers from scatter
        %use isoutlier function
scenario_2_days_to_plot_inhibited_tmp =scenario_2_plotdata.Days_to_plot;
scenario_2_plotdata_inhdivout_tmp = scenario_2_plotdata.InhibitedRadius_to_plot./scenario_2_plotdata.OuterRadius_to_plot;
[scenario_2_days_to_plot_inhdivout,scenario_2_plotdata_inhdivout] = function_removeoutliers_for_plot(scenario_2_days_to_plot_inhibited_tmp,scenario_2_plotdata_inhdivout_tmp);
 s2 = scatter(scenario_2_days_to_plot_inhdivout,scenario_2_plotdata_inhdivout,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Inhibited Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenarion_2_inhdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenarion_2_inhdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenarion_2_inhdivout_1'  '.png'])
print([filepath_save_figs  'Scenarion_2_inhdivout_1'],'-depsc2','-painters')

% BOX - nec/outer

figure
hold on
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 2),...
    scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2)./scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 4),...
    scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4)./scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 6),...
    scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6)./scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 8),...
    scenario_2_plotdata.NecroticRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8)./scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
        % exclude outliers from scatter
        %use isoutlier function
scenario_2_days_to_plot_necrotic_tmp =scenario_2_plotdata.Days_to_plot;
scenario_2_plotdata_necdivout_tmp = scenario_2_plotdata.NecroticRadius_to_plot./scenario_2_plotdata.OuterRadius_to_plot;
[scenario_2_days_to_plot_necdivout,scenario_2_plotdata_necdivout] = function_removeoutliers_for_plot(scenario_2_days_to_plot_necrotic_tmp,scenario_2_plotdata_necdivout_tmp);
 s2 = scatter(scenario_2_days_to_plot_necdivout,scenario_2_plotdata_necdivout,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Necrotic Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenarion_2_necdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenarion_2_necdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenarion_2_necdivout_1'  '.png'])
print([filepath_save_figs  'Scenarion_2_necdivout_1'],'-depsc2','-painters')


% BOX - PIM/outer

figure
hold on
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 2),...
    scenario_2_plotdata.PIMOuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2)./scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 2),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 4),...
    scenario_2_plotdata.PIMOuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4)./scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 4),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 6),...
    scenario_2_plotdata.PIMOuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6)./scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 6),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_2_plotdata.Days_to_plot(scenario_2_plotdata.Days_to_plot == 8),...
    scenario_2_plotdata.PIMOuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8)./scenario_2_plotdata.OuterRadius_to_plot(scenario_2_plotdata.Days_to_plot == 8),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])

        % exclude outliers from scatter
        %use isoutlier function
scenario_2_days_to_plot_pimouter_tmp =scenario_2_plotdata.Days_to_plot;
scenario_2_plotdata_pimouterdivout_tmp = scenario_2_plotdata.PIMOuterRadius_to_plot./scenario_2_plotdata.OuterRadius_to_plot;
[scenario_2_days_to_plot_pimouterdivout,scenario_2_plotdata_pimouterdivout] = function_removeoutliers_for_plot(scenario_2_days_to_plot_pimouter_tmp,scenario_2_plotdata_pimouterdivout_tmp);
 s2 = scatter(scenario_2_days_to_plot_pimouterdivout,scenario_2_plotdata_pimouterdivout,'filled','MarkerFaceColor',[0,1,1],'jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Hypoxic Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenarion_2_pimdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenarion_2_pimdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenarion_2_pimdivout_1'  '.png'])
print([filepath_save_figs  'Scenarion_2_pimdivout_1'],'-depsc2','-painters')



%% 3) PLOTS - Scenario 3 - Deoxygenation - 21% to Day 2 then 2%
scenario_3_conditions = [17,11,12,13,14];

[scenario_3_plotdata.OuterRadius_to_plot,...
    scenario_3_plotdata.InhibitedRadius_to_plot,...
    scenario_3_plotdata.NecroticRadius_to_plot,...
    scenario_3_plotdata.PIMOuterRadius_to_plot,...
    scenario_3_plotdata.PIMInnerRadius_to_plot,...
    scenario_3_plotdata.Days_to_plot,...
    scenario_3_plotdata.OuterRadius_mean_to_plot,...
    scenario_3_plotdata.InhibitedRadius_mean_to_plot,...
    scenario_3_plotdata.NecroticRadius_mean_to_plot,...
    scenario_3_plotdata.PIMOuterRadius_mean_to_plot,...
    scenario_3_plotdata.PIMInnerRadius_mean_to_plot,...
    scenario_3_plotdata.Days_unique_to_plot,...
    scenario_3_plotdata.OuterRadius_std_to_plot,...
    scenario_3_plotdata.InhibitedRadius_std_to_plot,...
    scenario_3_plotdata.NecroticRadius_std_to_plot,...
    scenario_3_plotdata.PIMOuterRadius_std_to_plot,...
    scenario_3_plotdata.PIMInnerRadius_std_to_plot,...
    scenario_3_plotdata.OuterRadius_count_to_plot,...
    scenario_3_plotdata.InhibitedRadius_count_to_plot,...
    scenario_3_plotdata.NecroticRadius_count_to_plot,...
    scenario_3_plotdata.PIMOuterRadius_count_to_plot,...
    scenario_3_plotdata.PIMInnerRadius_count_to_plot]  =  function_plotdata_scenario(scenario_3_conditions,...
    confocal_WM793b.AllData.Condition,...
    confocal_WM793b.AllData.OuterRadius,...
    confocal_WM793b.AllData.InhibitedRadius,...
    confocal_WM793b.AllData.NecroticRadius,...
    confocal_WM793b.AllData.PIMOuterRadius,...
    confocal_WM793b.AllData.PIMInnerRadius,...
    confocal_WM793b.AllData.Day);


% BOX PLOTS
figure
hold on
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 2),scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 3),scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 4),scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 6),scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 8),scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
hold on
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 2),scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 3),scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 4),scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 6),scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 8),scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
hold on
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 2),scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 3),scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 4),scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 6),scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 8),scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')

% exclude outliers from scatter
%use isoutlier function
scenario_3_days_to_plot_outer_tmp =scenario_3_plotdata.Days_to_plot;
scenario_3_plotdata_outer_tmp =  scenario_3_plotdata.OuterRadius_to_plot;
scenario_3_days_to_plot_inhibited_tmp =scenario_3_plotdata.Days_to_plot;
scenario_3_plotdata_inhibited_tmp = scenario_3_plotdata.InhibitedRadius_to_plot;
scenario_3_days_to_plot_necrotic_tmp = scenario_3_plotdata.Days_to_plot;
scenario_3_plotdata_necrotic_tmp = scenario_3_plotdata.NecroticRadius_to_plot;

[scenario_3_days_to_plot_outer,scenario_3_plotdata_outer] = function_removeoutliers_for_plot(scenario_3_days_to_plot_outer_tmp,scenario_3_plotdata_outer_tmp);
[scenario_3_days_to_plot_inhibited,scenario_3_plotdata_inhibited] = function_removeoutliers_for_plot(scenario_3_days_to_plot_inhibited_tmp,scenario_3_plotdata_inhibited_tmp);
[scenario_3_days_to_plot_necrotic,scenario_3_plotdata_necrotic] = function_removeoutliers_for_plot(scenario_3_days_to_plot_necrotic_tmp,scenario_3_plotdata_necrotic_tmp);

s3 = scatter(scenario_3_days_to_plot_outer,scenario_3_plotdata_outer,'g','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s2 = scatter(scenario_3_days_to_plot_inhibited,scenario_3_plotdata_inhibited,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(scenario_3_days_to_plot_necrotic,scenario_3_plotdata_necrotic,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);

alpha(s3,0.2)
alpha(s2,0.2)
alpha(s1,0.2)
xlim([0,8.5])
ylim([0,400])
xlabel('time [days]')
ylabel('Radius [microns]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_3_OIN_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_3_OIN_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_3_OIN_1'  '.png'])
print([filepath_save_figs  'Scenario_3_OIN_1'],'-depsc2','-painters')

% BOX PLOT WITH PIM

figure
hold on
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 2),scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 3),scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 4),scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 6),scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 8),scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
hold on
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 2),scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 3),scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 4),scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 6),scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 8),scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
hold on
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 2),scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 3),scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 4),scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 6),scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 8),scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
hold on
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 2),scenario_3_plotdata.PIMOuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 3),scenario_3_plotdata.PIMOuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 4),scenario_3_plotdata.PIMOuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 6),scenario_3_plotdata.PIMOuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 8),scenario_3_plotdata.PIMOuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])

% exclude outliers from scatter
%use isoutlier function
scenario_3_days_to_plot_outer_tmp =scenario_3_plotdata.Days_to_plot;
scenario_3_plotdata_outer_tmp =  scenario_3_plotdata.OuterRadius_to_plot;
scenario_3_days_to_plot_inhibited_tmp =scenario_3_plotdata.Days_to_plot;
scenario_3_plotdata_inhibited_tmp = scenario_3_plotdata.InhibitedRadius_to_plot;
scenario_3_days_to_plot_necrotic_tmp = scenario_3_plotdata.Days_to_plot;
scenario_3_plotdata_necrotic_tmp = scenario_3_plotdata.NecroticRadius_to_plot;
scenario_3_days_to_plot_PIMOuter_tmp = scenario_3_plotdata.Days_to_plot;
scenario_3_plotdata_PIMOuter_tmp = scenario_3_plotdata.PIMOuterRadius_to_plot;

[scenario_3_days_to_plot_outer,scenario_3_plotdata_outer] = function_removeoutliers_for_plot(scenario_3_days_to_plot_outer_tmp,scenario_3_plotdata_outer_tmp);
[scenario_3_days_to_plot_inhibited,scenario_3_plotdata_inhibited] = function_removeoutliers_for_plot(scenario_3_days_to_plot_inhibited_tmp,scenario_3_plotdata_inhibited_tmp);
[scenario_3_days_to_plot_necrotic,scenario_3_plotdata_necrotic] = function_removeoutliers_for_plot(scenario_3_days_to_plot_necrotic_tmp,scenario_3_plotdata_necrotic_tmp);
[scenario_3_days_to_plot_PIMOuter,scenario_3_plotdata_PIMOuter] = function_removeoutliers_for_plot(scenario_3_days_to_plot_PIMOuter_tmp,scenario_3_plotdata_PIMOuter_tmp);

s3 = scatter(scenario_3_days_to_plot_outer,scenario_3_plotdata_outer,'g','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s2 = scatter(scenario_3_days_to_plot_inhibited,scenario_3_plotdata_inhibited,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(scenario_3_days_to_plot_necrotic,scenario_3_plotdata_necrotic,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s0 = scatter(scenario_3_days_to_plot_PIMOuter,scenario_3_plotdata_PIMOuter,'filled','MarkerFaceColor',[0,1,1],'jitter', 'on', 'jitterAmount', jitter_plot);

alpha(s3,0.2)
alpha(s2,0.2)
alpha(s1,0.2)
alpha(s0,0.2)
xlim([0,8.5])
ylim([0,400])
xlabel('time [days]')
ylabel('Radius [microns]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_3_OINP_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_3_OINP_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_3_OINP_1'  '.png'])
print([filepath_save_figs  'Scenario_3_OINP_1'],'-depsc2','-painters')




% BOX PLOTS for inhibited radius/outer radius, necrotic radius/outer radius,
% pim outer radius/inhibited radius

% BOX - inh/outer

figure
hold on
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 2),...
    scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 3),...
    scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 4),...
    scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 6),...
    scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 8),...
    scenario_3_plotdata.InhibitedRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')

        % exclude outliers from scatter
        %use isoutlier function
scenario_3_days_to_plot_inhibited_tmp =scenario_3_plotdata.Days_to_plot;
scenario_3_plotdata_inhdivout_tmp = scenario_3_plotdata.InhibitedRadius_to_plot./scenario_3_plotdata.OuterRadius_to_plot;
[scenario_3_days_to_plot_inhdivout,scenario_3_plotdata_inhdivout] = function_removeoutliers_for_plot(scenario_3_days_to_plot_inhibited_tmp,scenario_3_plotdata_inhdivout_tmp);
 s2 = scatter(scenario_3_days_to_plot_inhdivout,scenario_3_plotdata_inhdivout,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Inhibited Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenarion_3_inhdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenarion_3_inhdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenarion_3_inhdivout_1'  '.png'])
print([filepath_save_figs  'Scenarion_3_inhdivout_1'],'-depsc2','-painters')

% BOX - nec/outer

figure
hold on
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 2),...
    scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 3),...
    scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 4),...
    scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 6),...
    scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 8),...
    scenario_3_plotdata.NecroticRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')

        % exclude outliers from scatter
        %use isoutlier function
scenario_3_days_to_plot_necrotic_tmp =scenario_3_plotdata.Days_to_plot;
scenario_3_plotdata_necdivout_tmp = scenario_3_plotdata.NecroticRadius_to_plot./scenario_3_plotdata.OuterRadius_to_plot;
[scenario_3_days_to_plot_necdivout,scenario_3_plotdata_necdivout] = function_removeoutliers_for_plot(scenario_3_days_to_plot_necrotic_tmp,scenario_3_plotdata_necdivout_tmp);
 s2 = scatter(scenario_3_days_to_plot_necdivout,scenario_3_plotdata_necdivout,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Necrotic Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenarion_3_necdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenarion_3_necdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenarion_3_necdivout_1'  '.png'])
print([filepath_save_figs  'Scenarion_3_necdivout_1'],'-depsc2','-painters')


% BOX - PIM/outer

figure
hold on
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 2),...
    scenario_3_plotdata.PIMOuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 2),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 3),...
     scenario_3_plotdata.PIMOuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 3),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 4),...
    scenario_3_plotdata.PIMOuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 4),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 6),...
    scenario_3_plotdata.PIMOuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 6),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_3_plotdata.Days_to_plot(scenario_3_plotdata.Days_to_plot == 8),...
    scenario_3_plotdata.PIMOuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8)./scenario_3_plotdata.OuterRadius_to_plot(scenario_3_plotdata.Days_to_plot == 8),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])

        % exclude outliers from scatter
        %use isoutlier function
scenario_3_days_to_plot_pimouter_tmp =scenario_3_plotdata.Days_to_plot;
scenario_3_plotdata_pimouterdivout_tmp = scenario_3_plotdata.PIMOuterRadius_to_plot./scenario_3_plotdata.OuterRadius_to_plot;
[scenario_3_days_to_plot_pimouterdivout,scenario_3_plotdata_pimouterdivout] = function_removeoutliers_for_plot(scenario_3_days_to_plot_pimouter_tmp,scenario_3_plotdata_pimouterdivout_tmp);
 s2 = scatter(scenario_3_days_to_plot_pimouterdivout,scenario_3_plotdata_pimouterdivout,'filled','MarkerFaceColor',[0,1,1],'jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Hypoxic Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenarion_3_pimdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenarion_3_pimdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenarion_3_pimdivout_1'  '.png'])
print([filepath_save_figs  'Scenarion_3_pimdivout_1'],'-depsc2','-painters')







%% 4) PLOTS - Scenario 4 - Re-oxygenation - 2% to Day 2 then 21%
scenario_4_conditions = [1,5,6,7,8];

[scenario_4_plotdata.OuterRadius_to_plot,...
    scenario_4_plotdata.InhibitedRadius_to_plot,...
    scenario_4_plotdata.NecroticRadius_to_plot,...
    scenario_4_plotdata.PIMOuterRadius_to_plot,...
    scenario_4_plotdata.PIMInnerRadius_to_plot,...
    scenario_4_plotdata.Days_to_plot,...
    scenario_4_plotdata.OuterRadius_mean_to_plot,...
    scenario_4_plotdata.InhibitedRadius_mean_to_plot,...
    scenario_4_plotdata.NecroticRadius_mean_to_plot,...
    scenario_4_plotdata.PIMOuterRadius_mean_to_plot,...
    scenario_4_plotdata.PIMInnerRadius_mean_to_plot,...
    scenario_4_plotdata.Days_unique_to_plot,...
    scenario_4_plotdata.OuterRadius_std_to_plot,...
    scenario_4_plotdata.InhibitedRadius_std_to_plot,...
    scenario_4_plotdata.NecroticRadius_std_to_plot,...
    scenario_4_plotdata.PIMOuterRadius_std_to_plot,...
    scenario_4_plotdata.PIMInnerRadius_std_to_plot,...
    scenario_4_plotdata.OuterRadius_count_to_plot,...
    scenario_4_plotdata.InhibitedRadius_count_to_plot,...
    scenario_4_plotdata.NecroticRadius_count_to_plot,...
    scenario_4_plotdata.PIMOuterRadius_count_to_plot,...
    scenario_4_plotdata.PIMInnerRadius_count_to_plot]  =  function_plotdata_scenario(scenario_4_conditions,...
    confocal_WM793b.AllData.Condition,...
    confocal_WM793b.AllData.OuterRadius,...
    confocal_WM793b.AllData.InhibitedRadius,...
    confocal_WM793b.AllData.NecroticRadius,...
    confocal_WM793b.AllData.PIMOuterRadius,...
    confocal_WM793b.AllData.PIMInnerRadius,...
    confocal_WM793b.AllData.Day);


% BOX PLOTS
figure
hold on
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 2),scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 3),scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 4),scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 6),scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 8),scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
hold on
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 2),scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 3),scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 4),scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 6),scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 8),scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
hold on
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 2),scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 3),scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 4),scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 6),scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 8),scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')

% exclude outliers from scatter
%use isoutlier function
scenario_4_days_to_plot_outer_tmp =scenario_4_plotdata.Days_to_plot;
scenario_4_plotdata_outer_tmp =  scenario_4_plotdata.OuterRadius_to_plot;
scenario_4_days_to_plot_inhibited_tmp =scenario_4_plotdata.Days_to_plot;
scenario_4_plotdata_inhibited_tmp = scenario_4_plotdata.InhibitedRadius_to_plot;
scenario_4_days_to_plot_necrotic_tmp = scenario_4_plotdata.Days_to_plot;
scenario_4_plotdata_necrotic_tmp = scenario_4_plotdata.NecroticRadius_to_plot;

[scenario_4_days_to_plot_outer,scenario_4_plotdata_outer] = function_removeoutliers_for_plot(scenario_4_days_to_plot_outer_tmp,scenario_4_plotdata_outer_tmp);
[scenario_4_days_to_plot_inhibited,scenario_4_plotdata_inhibited] = function_removeoutliers_for_plot(scenario_4_days_to_plot_inhibited_tmp,scenario_4_plotdata_inhibited_tmp);
[scenario_4_days_to_plot_necrotic,scenario_4_plotdata_necrotic] = function_removeoutliers_for_plot(scenario_4_days_to_plot_necrotic_tmp,scenario_4_plotdata_necrotic_tmp);

s3 = scatter(scenario_4_days_to_plot_outer,scenario_4_plotdata_outer,'g','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s2 = scatter(scenario_4_days_to_plot_inhibited,scenario_4_plotdata_inhibited,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(scenario_4_days_to_plot_necrotic,scenario_4_plotdata_necrotic,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);

alpha(s3,0.2)
alpha(s2,0.2)
alpha(s1,0.2)
xlim([0,8.5])
ylim([0,400])
xlabel('time [days]')
ylabel('Radius [microns]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_4_OIN_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_4_OIN_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_4_OIN_1'  '.png'])
print([filepath_save_figs  'Scenario_4_OIN_1'],'-depsc2','-painters')

% BOX PLOT WITH PIM

figure
hold on
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 2),scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 3),scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 4),scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 6),scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 8),scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
hold on
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 2),scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 3),scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 4),scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 6),scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 8),scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
hold on
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 2),scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 3),scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 4),scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 6),scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 8),scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
hold on
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 2),scenario_4_plotdata.PIMOuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 3),scenario_4_plotdata.PIMOuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 4),scenario_4_plotdata.PIMOuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 6),scenario_4_plotdata.PIMOuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 8),scenario_4_plotdata.PIMOuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])

% exclude outliers from scatter
%use isoutlier function
scenario_4_days_to_plot_outer_tmp =scenario_4_plotdata.Days_to_plot;
scenario_4_plotdata_outer_tmp =  scenario_4_plotdata.OuterRadius_to_plot;
scenario_4_days_to_plot_inhibited_tmp =scenario_4_plotdata.Days_to_plot;
scenario_4_plotdata_inhibited_tmp = scenario_4_plotdata.InhibitedRadius_to_plot;
scenario_4_days_to_plot_necrotic_tmp = scenario_4_plotdata.Days_to_plot;
scenario_4_plotdata_necrotic_tmp = scenario_4_plotdata.NecroticRadius_to_plot;
scenario_4_days_to_plot_PIMOuter_tmp = scenario_4_plotdata.Days_to_plot;
scenario_4_plotdata_PIMOuter_tmp = scenario_4_plotdata.PIMOuterRadius_to_plot;

[scenario_4_days_to_plot_outer,scenario_4_plotdata_outer] = function_removeoutliers_for_plot(scenario_4_days_to_plot_outer_tmp,scenario_4_plotdata_outer_tmp);
[scenario_4_days_to_plot_inhibited,scenario_4_plotdata_inhibited] = function_removeoutliers_for_plot(scenario_4_days_to_plot_inhibited_tmp,scenario_4_plotdata_inhibited_tmp);
[scenario_4_days_to_plot_necrotic,scenario_4_plotdata_necrotic] = function_removeoutliers_for_plot(scenario_4_days_to_plot_necrotic_tmp,scenario_4_plotdata_necrotic_tmp);
[scenario_4_days_to_plot_PIMOuter,scenario_4_plotdata_PIMOuter] = function_removeoutliers_for_plot(scenario_4_days_to_plot_PIMOuter_tmp,scenario_4_plotdata_PIMOuter_tmp);

s3 = scatter(scenario_4_days_to_plot_outer,scenario_4_plotdata_outer,'g','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s2 = scatter(scenario_4_days_to_plot_inhibited,scenario_4_plotdata_inhibited,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(scenario_4_days_to_plot_necrotic,scenario_4_plotdata_necrotic,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s0 = scatter(scenario_4_days_to_plot_PIMOuter,scenario_4_plotdata_PIMOuter,'filled','MarkerFaceColor',[0,1,1],'jitter', 'on', 'jitterAmount', jitter_plot);

alpha(s3,0.2)
alpha(s2,0.2)
alpha(s1,0.2)
alpha(s0,0.2)
xlim([0,8.5])
ylim([0,400])
xlabel('time [days]')
ylabel('Radius [microns]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_4_OINP_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_4_OINP_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_4_OINP_1'  '.png'])
print([filepath_save_figs  'Scenario_4_OINP_1'],'-depsc2','-painters')




% BOX PLOTS for inhibited radius/outer radius, necrotic radius/outer radius,
% pim outer radius/inhibited radius

% BOX - inh/outer

figure
hold on
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 2),...
    scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 3),...
    scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 4),...
    scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 6),...
    scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 8),...
    scenario_4_plotdata.InhibitedRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')

        % exclude outliers from scatter
        %use isoutlier function
scenario_4_days_to_plot_inhibited_tmp =scenario_4_plotdata.Days_to_plot;
scenario_4_plotdata_inhdivout_tmp = scenario_4_plotdata.InhibitedRadius_to_plot./scenario_4_plotdata.OuterRadius_to_plot;
[scenario_4_days_to_plot_inhdivout,scenario_4_plotdata_inhdivout] = function_removeoutliers_for_plot(scenario_4_days_to_plot_inhibited_tmp,scenario_4_plotdata_inhdivout_tmp);
 s2 = scatter(scenario_4_days_to_plot_inhdivout,scenario_4_plotdata_inhdivout,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Inhibited Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenarion_4_inhdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenarion_4_inhdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenarion_4_inhdivout_1'  '.png'])
print([filepath_save_figs  'Scenarion_4_inhdivout_1'],'-depsc2','-painters')

% BOX - nec/outer

figure
hold on
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 2),...
    scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 3),...
    scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 4),...
    scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 6),...
    scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 8),...
    scenario_4_plotdata.NecroticRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')

        % exclude outliers from scatter
        %use isoutlier function
scenario_4_days_to_plot_necrotic_tmp =scenario_4_plotdata.Days_to_plot;
scenario_4_plotdata_necdivout_tmp = scenario_4_plotdata.NecroticRadius_to_plot./scenario_4_plotdata.OuterRadius_to_plot;
[scenario_4_days_to_plot_necdivout,scenario_4_plotdata_necdivout] = function_removeoutliers_for_plot(scenario_4_days_to_plot_necrotic_tmp,scenario_4_plotdata_necdivout_tmp);
 s2 = scatter(scenario_4_days_to_plot_necdivout,scenario_4_plotdata_necdivout,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Necrotic Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenarion_4_necdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenarion_4_necdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenarion_4_necdivout_1'  '.png'])
print([filepath_save_figs  'Scenarion_4_necdivout_1'],'-depsc2','-painters')


% BOX - PIM/outer

figure
hold on
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 2),...
    scenario_4_plotdata.PIMOuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 2),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 3),...
     scenario_4_plotdata.PIMOuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 3),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 4),...
    scenario_4_plotdata.PIMOuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 4),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 6),...
    scenario_4_plotdata.PIMOuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 6),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_4_plotdata.Days_to_plot(scenario_4_plotdata.Days_to_plot == 8),...
    scenario_4_plotdata.PIMOuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8)./scenario_4_plotdata.OuterRadius_to_plot(scenario_4_plotdata.Days_to_plot == 8),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])

        % exclude outliers from scatter
        %use isoutlier function
scenario_4_days_to_plot_pimouter_tmp =scenario_4_plotdata.Days_to_plot;
scenario_4_plotdata_pimouterdivout_tmp = scenario_4_plotdata.PIMOuterRadius_to_plot./scenario_4_plotdata.OuterRadius_to_plot;
[scenario_4_days_to_plot_pimouterdivout,scenario_4_plotdata_pimouterdivout] = function_removeoutliers_for_plot(scenario_4_days_to_plot_pimouter_tmp,scenario_4_plotdata_pimouterdivout_tmp);
 s2 = scatter(scenario_4_days_to_plot_pimouterdivout,scenario_4_plotdata_pimouterdivout,'filled','MarkerFaceColor',[0,1,1],'jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Hypoxic Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenarion_4_pimdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenarion_4_pimdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenarion_4_pimdivout_1'  '.png'])
print([filepath_save_figs  'Scenarion_4_pimdivout_1'],'-depsc2','-painters')

%% 5) PLOTS - Scenario 5 - Re-oxygenation - 2% to Day 4 then 21%
scenario_5_conditions = [1,2,9,10];

[scenario_5_plotdata.OuterRadius_to_plot,...
    scenario_5_plotdata.InhibitedRadius_to_plot,...
    scenario_5_plotdata.NecroticRadius_to_plot,...
    scenario_5_plotdata.PIMOuterRadius_to_plot,...
    scenario_5_plotdata.PIMInnerRadius_to_plot,...
    scenario_5_plotdata.Days_to_plot,...
    scenario_5_plotdata.OuterRadius_mean_to_plot,...
    scenario_5_plotdata.InhibitedRadius_mean_to_plot,...
    scenario_5_plotdata.NecroticRadius_mean_to_plot,...
    scenario_5_plotdata.PIMOuterRadius_mean_to_plot,...
    scenario_5_plotdata.PIMInnerRadius_mean_to_plot,...
    scenario_5_plotdata.Days_unique_to_plot,...
    scenario_5_plotdata.OuterRadius_std_to_plot,...
    scenario_5_plotdata.InhibitedRadius_std_to_plot,...
    scenario_5_plotdata.NecroticRadius_std_to_plot,...
    scenario_5_plotdata.PIMOuterRadius_std_to_plot,...
    scenario_5_plotdata.PIMInnerRadius_std_to_plot,...
    scenario_5_plotdata.OuterRadius_count_to_plot,...
    scenario_5_plotdata.InhibitedRadius_count_to_plot,...
    scenario_5_plotdata.NecroticRadius_count_to_plot,...
    scenario_5_plotdata.PIMOuterRadius_count_to_plot,...
    scenario_5_plotdata.PIMInnerRadius_count_to_plot]  =  function_plotdata_scenario(scenario_5_conditions,...
    confocal_WM793b.AllData.Condition,...
    confocal_WM793b.AllData.OuterRadius,...
    confocal_WM793b.AllData.InhibitedRadius,...
    confocal_WM793b.AllData.NecroticRadius,...
    confocal_WM793b.AllData.PIMOuterRadius,...
    confocal_WM793b.AllData.PIMInnerRadius,...
    confocal_WM793b.AllData.Day);


% BOX PLOTS
figure
hold on
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 2),scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 3),scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 4),scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 6),scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 8),scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
hold on
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 2),scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 3),scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 4),scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 6),scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 8),scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
hold on
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 2),scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 3),scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 4),scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 6),scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 8),scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')

% exclude outliers from scatter
%use isoutlier function
scenario_5_days_to_plot_outer_tmp =scenario_5_plotdata.Days_to_plot;
scenario_5_plotdata_outer_tmp =  scenario_5_plotdata.OuterRadius_to_plot;
scenario_5_days_to_plot_inhibited_tmp =scenario_5_plotdata.Days_to_plot;
scenario_5_plotdata_inhibited_tmp = scenario_5_plotdata.InhibitedRadius_to_plot;
scenario_5_days_to_plot_necrotic_tmp = scenario_5_plotdata.Days_to_plot;
scenario_5_plotdata_necrotic_tmp = scenario_5_plotdata.NecroticRadius_to_plot;

[scenario_5_days_to_plot_outer,scenario_5_plotdata_outer] = function_removeoutliers_for_plot(scenario_5_days_to_plot_outer_tmp,scenario_5_plotdata_outer_tmp);
[scenario_5_days_to_plot_inhibited,scenario_5_plotdata_inhibited] = function_removeoutliers_for_plot(scenario_5_days_to_plot_inhibited_tmp,scenario_5_plotdata_inhibited_tmp);
[scenario_5_days_to_plot_necrotic,scenario_5_plotdata_necrotic] = function_removeoutliers_for_plot(scenario_5_days_to_plot_necrotic_tmp,scenario_5_plotdata_necrotic_tmp);

s3 = scatter(scenario_5_days_to_plot_outer,scenario_5_plotdata_outer,'g','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s2 = scatter(scenario_5_days_to_plot_inhibited,scenario_5_plotdata_inhibited,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(scenario_5_days_to_plot_necrotic,scenario_5_plotdata_necrotic,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);

alpha(s3,0.2)
alpha(s2,0.2)
alpha(s1,0.2)
xlim([0,8.5])
ylim([0,400])
xlabel('time [days]')
ylabel('Radius [microns]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_5_OIN_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_5_OIN_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_5_OIN_1'  '.png'])
print([filepath_save_figs  'Scenario_5_OIN_1'],'-depsc2','-painters')

% BOX PLOT WITH PIM

figure
hold on
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 2),scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 3),scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 4),scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 6),scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 8),scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8),'BoxFaceColor','g','JitterOutliers', 'on','MarkerColor','g')
hold on
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 2),scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 3),scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 4),scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 6),scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 8),scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
hold on
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 2),scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 3),scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 4),scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 6),scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 8),scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
hold on
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 2),scenario_5_plotdata.PIMOuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 3),scenario_5_plotdata.PIMOuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 4),scenario_5_plotdata.PIMOuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 6),scenario_5_plotdata.PIMOuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 8),scenario_5_plotdata.PIMOuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])

% exclude outliers from scatter
%use isoutlier function
scenario_5_days_to_plot_outer_tmp =scenario_5_plotdata.Days_to_plot;
scenario_5_plotdata_outer_tmp =  scenario_5_plotdata.OuterRadius_to_plot;
scenario_5_days_to_plot_inhibited_tmp =scenario_5_plotdata.Days_to_plot;
scenario_5_plotdata_inhibited_tmp = scenario_5_plotdata.InhibitedRadius_to_plot;
scenario_5_days_to_plot_necrotic_tmp = scenario_5_plotdata.Days_to_plot;
scenario_5_plotdata_necrotic_tmp = scenario_5_plotdata.NecroticRadius_to_plot;
scenario_5_days_to_plot_PIMOuter_tmp = scenario_5_plotdata.Days_to_plot;
scenario_5_plotdata_PIMOuter_tmp = scenario_5_plotdata.PIMOuterRadius_to_plot;

[scenario_5_days_to_plot_outer,scenario_5_plotdata_outer] = function_removeoutliers_for_plot(scenario_5_days_to_plot_outer_tmp,scenario_5_plotdata_outer_tmp);
[scenario_5_days_to_plot_inhibited,scenario_5_plotdata_inhibited] = function_removeoutliers_for_plot(scenario_5_days_to_plot_inhibited_tmp,scenario_5_plotdata_inhibited_tmp);
[scenario_5_days_to_plot_necrotic,scenario_5_plotdata_necrotic] = function_removeoutliers_for_plot(scenario_5_days_to_plot_necrotic_tmp,scenario_5_plotdata_necrotic_tmp);
[scenario_5_days_to_plot_PIMOuter,scenario_5_plotdata_PIMOuter] = function_removeoutliers_for_plot(scenario_5_days_to_plot_PIMOuter_tmp,scenario_5_plotdata_PIMOuter_tmp);

s3 = scatter(scenario_5_days_to_plot_outer,scenario_5_plotdata_outer,'g','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s2 = scatter(scenario_5_days_to_plot_inhibited,scenario_5_plotdata_inhibited,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(scenario_5_days_to_plot_necrotic,scenario_5_plotdata_necrotic,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s0 = scatter(scenario_5_days_to_plot_PIMOuter,scenario_5_plotdata_PIMOuter,'filled','MarkerFaceColor',[0,1,1],'jitter', 'on', 'jitterAmount', jitter_plot);

alpha(s3,0.2)
alpha(s2,0.2)
alpha(s1,0.2)
alpha(s0,0.2)
xlim([0,8.5])
ylim([0,400])
xlabel('time [days]')
ylabel('Radius [microns]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenario_5_OINP_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_5_OINP_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_5_OINP_1'  '.png'])
print([filepath_save_figs  'Scenario_5_OINP_1'],'-depsc2','-painters')


% BOX PLOTS for inhibited radius/outer radius, necrotic radius/outer radius,
% pim outer radius/inhibited radius

% BOX - inh/outer

figure
hold on
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 2),...
    scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 3),...
    scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 4),...
    scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 6),...
    scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 8),...
    scenario_5_plotdata.InhibitedRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8),'BoxFaceColor','m','JitterOutliers', 'on','MarkerColor','m')

        % exclude outliers from scatter
        %use isoutlier function
scenario_5_days_to_plot_inhibited_tmp =scenario_5_plotdata.Days_to_plot;
scenario_5_plotdata_inhdivout_tmp = scenario_5_plotdata.InhibitedRadius_to_plot./scenario_5_plotdata.OuterRadius_to_plot;
[scenario_5_days_to_plot_inhdivout,scenario_5_plotdata_inhdivout] = function_removeoutliers_for_plot(scenario_5_days_to_plot_inhibited_tmp,scenario_5_plotdata_inhdivout_tmp);
 s2 = scatter(scenario_5_days_to_plot_inhdivout,scenario_5_plotdata_inhdivout,'m','filled','jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Inhibited Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenarion_5_inhdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenarion_5_inhdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenarion_5_inhdivout_1'  '.png'])
print([filepath_save_figs  'Scenarion_5_inhdivout_1'],'-depsc2','-painters')

% BOX - nec/outer

figure
hold on
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 2),...
    scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 3),...
    scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 4),...
    scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 6),...
    scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 8),...
    scenario_5_plotdata.NecroticRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8),'BoxFaceColor','k','JitterOutliers', 'on','MarkerColor','k')

        % exclude outliers from scatter
        %use isoutlier function
scenario_5_days_to_plot_necrotic_tmp =scenario_5_plotdata.Days_to_plot;
scenario_5_plotdata_necdivout_tmp = scenario_5_plotdata.NecroticRadius_to_plot./scenario_5_plotdata.OuterRadius_to_plot;
[scenario_5_days_to_plot_necdivout,scenario_5_plotdata_necdivout] = function_removeoutliers_for_plot(scenario_5_days_to_plot_necrotic_tmp,scenario_5_plotdata_necdivout_tmp);
 s2 = scatter(scenario_5_days_to_plot_necdivout,scenario_5_plotdata_necdivout,'k','filled','jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Necrotic Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenarion_5_necdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenarion_5_necdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenarion_5_necdivout_1'  '.png'])
print([filepath_save_figs  'Scenarion_5_necdivout_1'],'-depsc2','-painters')


% BOX - PIM/outer

figure
hold on
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 2),...
    scenario_5_plotdata.PIMOuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 2),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 3),...
     scenario_5_plotdata.PIMOuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 3),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 4),...
    scenario_5_plotdata.PIMOuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 4),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 6),...
    scenario_5_plotdata.PIMOuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 6),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])
boxchart(scenario_5_plotdata.Days_to_plot(scenario_5_plotdata.Days_to_plot == 8),...
    scenario_5_plotdata.PIMOuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8)./scenario_5_plotdata.OuterRadius_to_plot(scenario_5_plotdata.Days_to_plot == 8),'BoxFaceColor',[0,1,1],'JitterOutliers', 'on','MarkerColor',[0,1,1])

        % exclude outliers from scatter
        %use isoutlier function
scenario_5_days_to_plot_pimouter_tmp =scenario_5_plotdata.Days_to_plot;
scenario_5_plotdata_pimouterdivout_tmp = scenario_5_plotdata.PIMOuterRadius_to_plot./scenario_5_plotdata.OuterRadius_to_plot;
[scenario_5_days_to_plot_pimouterdivout,scenario_5_plotdata_pimouterdivout] = function_removeoutliers_for_plot(scenario_5_days_to_plot_pimouter_tmp,scenario_5_plotdata_pimouterdivout_tmp);
 s2 = scatter(scenario_5_days_to_plot_pimouterdivout,scenario_5_plotdata_pimouterdivout,'filled','MarkerFaceColor',[0,1,1],'jitter', 'on', 'jitterAmount', jitter_plot);
 alpha(s2,0.2)
xlim([0,8.5])
ylim([0,1])
xlabel('time [days]')
ylabel('Hypoxic Radius/ Outer radius [-]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'Scenarion_5_pimdivout_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenarion_5_pimdivout_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenarion_5_pimdivout_1'  '.png'])
print([filepath_save_figs  'Scenarion_5_pimdivout_1'],'-depsc2','-painters')



%% 6) PLOTS - Compare Scenario 1 and 2


figure
hold on
scatter(scenario_1_plotdata.OuterRadius_to_plot,...
    scenario_1_plotdata.InhibitedRadius_to_plot./scenario_1_plotdata.OuterRadius_to_plot,'filled','MarkerFaceColor','r','MarkerEdgeColor','r');
scatter(scenario_2_plotdata.OuterRadius_to_plot,...
    scenario_2_plotdata.InhibitedRadius_to_plot./scenario_2_plotdata.OuterRadius_to_plot,'filled','MarkerFaceColor','b','MarkerEdgeColor','b');
xlim([0,400])
ylim([0,1])
legend('Normoxia','Hypoxia','Location','southwest')
title('Inhibited radius / Outer radius')
xlabel('Outer radius')
ylabel('Inhibited radius / Outer radius')
box on
grid on
saveas(gcf, [filepath_save_figs  'Scenario_Comparison_1_2_inhdivout'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_Comparison_1_2_inhdivout'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_Comparison_1_2_inhdivout'  '.png'])
print([filepath_save_figs  'Scenario_Comparison_1_2_inhdivout'],'-depsc2','-painters')

figure
hold on
scatter(scenario_1_plotdata.OuterRadius_to_plot,...
    scenario_1_plotdata.NecroticRadius_to_plot./scenario_1_plotdata.OuterRadius_to_plot,'filled','MarkerFaceColor','r','MarkerEdgeColor','r');
scatter(scenario_2_plotdata.OuterRadius_to_plot,...
    scenario_2_plotdata.NecroticRadius_to_plot./scenario_2_plotdata.OuterRadius_to_plot,'filled','MarkerFaceColor','b','MarkerEdgeColor','b');
xlim([0,400])
ylim([0,1])
legend('Normoxia','Hypoxia','Location','southwest')
title('Necrotic radius / Outer radius')
xlabel('Outer radius')
ylabel('Necrotic radius / Outer radius')
box on
grid on
saveas(gcf, [filepath_save_figs  'Scenario_Comparison_1_2_necdivout'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_Comparison_1_2_necdivout'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_Comparison_1_2_necdivout'  '.png'])
print([filepath_save_figs  'Scenario_Comparison_1_2_necdivout'],'-depsc2','-painters')


figure
hold on
scatter(scenario_1_plotdata.OuterRadius_to_plot,...
    scenario_1_plotdata.PIMOuterRadius_to_plot./scenario_1_plotdata.OuterRadius_to_plot,'filled','MarkerFaceColor','r','MarkerEdgeColor','r');
scatter(scenario_2_plotdata.OuterRadius_to_plot,...
    scenario_2_plotdata.PIMOuterRadius_to_plot./scenario_2_plotdata.OuterRadius_to_plot,'filled','MarkerFaceColor','b','MarkerEdgeColor','b');
xlim([0,400])
ylim([0,1])
legend('Normoxia','Hypoxia','Location','southwest')
title('Hypoxic radius / Outer radius')
xlabel('Outer radius')
ylabel('Hypoxic radius / Outer radius')
box on
grid on

saveas(gcf, [filepath_save_figs  'Scenario_Comparison_1_2_pimdivout'  '.fig'])
saveas(gcf, [filepath_save_figs  'Scenario_Comparison_1_2_pimdivout'  '.pdf'])
saveas(gcf, [filepath_save_figs  'Scenario_Comparison_1_2_pimdivout'  '.png'])
print([filepath_save_figs  'Scenario_Comparison_1_2_pimdivout'],'-depsc2','-painters')



%% 7) Save data for profile likelihood analysis
close all


%% 7.1) Save outer radius data

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.OuterRadius;
Outlier = zeros(length(confocal_WM793b.AllData.OuterRadius),1);

expected_number_of_outliers = length(Radius) - length(scenario_1_plotdata_outer) - length(scenario_2_plotdata_outer);

loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_1_conditions)>0
        % if in list keep
        Scenario(i) = 1;
        if sum((Day(i)== scenario_1_days_to_plot_outer).*(Radius(i)==scenario_1_plotdata_outer)) > 0
            Outlier(i) = 0;
        else
            Outlier(i) = 1;
        end
    elseif  sum(Condition(i) == scenario_2_conditions)>0
        Scenario(i) = 2;
        if sum((Day(i)== scenario_2_days_to_plot_outer).*(Radius(i)==scenario_2_plotdata_outer)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    elseif  sum(Condition(i) == scenario_3_conditions)>0
        Scenario(i) = 3;
        if sum((Day(i)== scenario_3_days_to_plot_outer).*(Radius(i)==scenario_3_plotdata_outer)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    elseif  sum(Condition(i) == scenario_4_conditions)>0
        Scenario(i) = 4;
        if sum((Day(i)== scenario_4_days_to_plot_outer).*(Radius(i)==scenario_4_plotdata_outer)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    elseif  sum(Condition(i) == scenario_5_conditions)>0
        Scenario(i) = 5;
        if sum((Day(i)== scenario_5_days_to_plot_outer).*(Radius(i)==scenario_5_plotdata_outer)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end
% sum(Outlier)


save([filepath_save 'WM793b_confocal_outer' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;

fname = fname_tmpscen(Scenario==1);
ExpID = ExpID_tmpscen(Scenario==1);
CellLine = CellLine_tmpscen(Scenario==1); 
Condition = Condition_tmpscen(Scenario==1);
Day = Day_tmpscen(Scenario==1);
Replicate = Replicate_tmpscen(Scenario==1);
Radius = Radius_tmpscen(Scenario==1);
Outlier = Outlier_tmpscen(Scenario==1);

save([filepath_save 'WM793b_confocal_scenario_1_outer' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');

fname = fname_tmpscen(Scenario==2);
ExpID = ExpID_tmpscen(Scenario==2);
CellLine = CellLine_tmpscen(Scenario==2); 
Condition = Condition_tmpscen(Scenario==2);
Day = Day_tmpscen(Scenario==2);
Replicate = Replicate_tmpscen(Scenario==2);
Radius = Radius_tmpscen(Scenario==2);
Outlier = Outlier_tmpscen(Scenario==2);

save([filepath_save 'WM793b_confocal_scenario_2_outer' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');

fname = fname_tmpscen(Scenario==3);
ExpID = ExpID_tmpscen(Scenario==3);
CellLine = CellLine_tmpscen(Scenario==3); 
Condition = Condition_tmpscen(Scenario==3);
Day = Day_tmpscen(Scenario==3);
Replicate = Replicate_tmpscen(Scenario==3);
Radius = Radius_tmpscen(Scenario==3);
Outlier = Outlier_tmpscen(Scenario==3);

save([filepath_save 'WM793b_confocal_scenario_3_outer' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


fname = fname_tmpscen(Scenario==4);
ExpID = ExpID_tmpscen(Scenario==4);
CellLine = CellLine_tmpscen(Scenario==4); 
Condition = Condition_tmpscen(Scenario==4);
Day = Day_tmpscen(Scenario==4);
Replicate = Replicate_tmpscen(Scenario==4);
Radius = Radius_tmpscen(Scenario==4);
Outlier = Outlier_tmpscen(Scenario==4);

save([filepath_save 'WM793b_confocal_scenario_4_outer' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');

fname = fname_tmpscen(Scenario==5);
ExpID = ExpID_tmpscen(Scenario==5);
CellLine = CellLine_tmpscen(Scenario==5); 
Condition = Condition_tmpscen(Scenario==5);
Day = Day_tmpscen(Scenario==5);
Replicate = Replicate_tmpscen(Scenario==5);
Radius = Radius_tmpscen(Scenario==5);
Outlier = Outlier_tmpscen(Scenario==5);

save([filepath_save 'WM793b_confocal_scenario_5_outer' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


%% 7.2) Save inhibited radius data

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.InhibitedRadius;
Outlier = zeros(length(confocal_WM793b.AllData.InhibitedRadius),1);

expected_number_of_outliers = length(Radius) - length(scenario_1_plotdata_inhibited) - length(scenario_2_plotdata_inhibited);

loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_1_conditions)>0
        Scenario(i) = 1;
        % if in list keep
        if sum((Day(i)== scenario_1_days_to_plot_inhibited).*(Radius(i)==scenario_1_plotdata_inhibited)) > 0
            Outlier(i) = 0;
        else
            Outlier(i) = 1;
        end
    elseif  sum(Condition(i) == scenario_2_conditions)>0
        Scenario(i) = 2;
        if sum((Day(i)== scenario_2_days_to_plot_inhibited).*(Radius(i)==scenario_2_plotdata_inhibited)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    elseif  sum(Condition(i) == scenario_3_conditions)>0
        Scenario(i) = 3;
        if sum((Day(i)== scenario_3_days_to_plot_inhibited).*(Radius(i)==scenario_3_plotdata_inhibited)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    elseif  sum(Condition(i) == scenario_4_conditions)>0
        Scenario(i) = 4;
        if sum((Day(i)== scenario_4_days_to_plot_inhibited).*(Radius(i)==scenario_4_plotdata_inhibited)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    elseif  sum(Condition(i) == scenario_5_conditions)>0
        Scenario(i) = 5;
        if sum((Day(i)== scenario_5_days_to_plot_inhibited).*(Radius(i)==scenario_5_plotdata_inhibited)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end
% sum(Outlier)


save([filepath_save 'WM793b_confocal_inhibited' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');

fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;


fname = fname_tmpscen(Scenario==1);
ExpID = ExpID_tmpscen(Scenario==1);
CellLine = CellLine_tmpscen(Scenario==1); 
Condition = Condition_tmpscen(Scenario==1);
Day = Day_tmpscen(Scenario==1);
Replicate = Replicate_tmpscen(Scenario==1);
Radius = Radius_tmpscen(Scenario==1);
Outlier = Outlier_tmpscen(Scenario==1);

save([filepath_save 'WM793b_confocal_scenario_1_inhibited' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');

fname = fname_tmpscen(Scenario==2);
ExpID = ExpID_tmpscen(Scenario==2);
CellLine = CellLine_tmpscen(Scenario==2); 
Condition = Condition_tmpscen(Scenario==2);
Day = Day_tmpscen(Scenario==2);
Replicate = Replicate_tmpscen(Scenario==2);
Radius = Radius_tmpscen(Scenario==2);
Outlier = Outlier_tmpscen(Scenario==2);

save([filepath_save 'WM793b_confocal_scenario_2_inhibited' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');

fname = fname_tmpscen(Scenario==3);
ExpID = ExpID_tmpscen(Scenario==3);
CellLine = CellLine_tmpscen(Scenario==3); 
Condition = Condition_tmpscen(Scenario==3);
Day = Day_tmpscen(Scenario==3);
Replicate = Replicate_tmpscen(Scenario==3);
Radius = Radius_tmpscen(Scenario==3);
Outlier = Outlier_tmpscen(Scenario==3);

save([filepath_save 'WM793b_confocal_scenario_3_inhibited' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


fname = fname_tmpscen(Scenario==4);
ExpID = ExpID_tmpscen(Scenario==4);
CellLine = CellLine_tmpscen(Scenario==4); 
Condition = Condition_tmpscen(Scenario==4);
Day = Day_tmpscen(Scenario==4);
Replicate = Replicate_tmpscen(Scenario==4);
Radius = Radius_tmpscen(Scenario==4);
Outlier = Outlier_tmpscen(Scenario==4);

save([filepath_save 'WM793b_confocal_scenario_4_inhibited' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


fname = fname_tmpscen(Scenario==5);
ExpID = ExpID_tmpscen(Scenario==5);
CellLine = CellLine_tmpscen(Scenario==5); 
Condition = Condition_tmpscen(Scenario==5);
Day = Day_tmpscen(Scenario==5);
Replicate = Replicate_tmpscen(Scenario==5);
Radius = Radius_tmpscen(Scenario==5);
Outlier = Outlier_tmpscen(Scenario==5);

save([filepath_save 'WM793b_confocal_scenario_5_inhibited' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


%% 7.3) Save necrotic radius data


fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.NecroticRadius;
Outlier = zeros(length(confocal_WM793b.AllData.NecroticRadius),1);

expected_number_of_outliers = length(Radius) - length(scenario_1_plotdata_necrotic) - length(scenario_2_plotdata_necrotic);

loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_1_conditions)>0
        Scenario(i) = 1;
        % if in list keep
        if sum((Day(i)== scenario_1_days_to_plot_necrotic).*(Radius(i)==scenario_1_plotdata_necrotic)) > 0
            Outlier(i) = 0;
        else
            Outlier(i) = 1;
        end
    elseif  sum(Condition(i) == scenario_2_conditions)>0
        Scenario(i) = 2;
        if sum((Day(i)== scenario_2_days_to_plot_necrotic).*(Radius(i)==scenario_2_plotdata_necrotic)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    elseif  sum(Condition(i) == scenario_3_conditions)>0
       Scenario(i) = 3; 
        if sum((Day(i)== scenario_3_days_to_plot_necrotic).*(Radius(i)==scenario_3_plotdata_necrotic)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    elseif  sum(Condition(i) == scenario_4_conditions)>0
        Scenario(i) = 4;
        if sum((Day(i)== scenario_4_days_to_plot_necrotic).*(Radius(i)==scenario_4_plotdata_necrotic)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    elseif  sum(Condition(i) == scenario_5_conditions)>0
        Scenario(i) = 5;
        if sum((Day(i)== scenario_5_days_to_plot_necrotic).*(Radius(i)==scenario_5_plotdata_necrotic)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end

% sum(Outlier)


save([filepath_save 'WM793b_confocal_necrotic' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');

fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;



fname = fname_tmpscen(Scenario==1);
ExpID = ExpID_tmpscen(Scenario==1);
CellLine = CellLine_tmpscen(Scenario==1); 
Condition = Condition_tmpscen(Scenario==1);
Day = Day_tmpscen(Scenario==1);
Replicate = Replicate_tmpscen(Scenario==1);
Radius = Radius_tmpscen(Scenario==1);
Outlier = Outlier_tmpscen(Scenario==1);

save([filepath_save 'WM793b_confocal_scenario_1_necrotic' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');

fname = fname_tmpscen(Scenario==2);
ExpID = ExpID_tmpscen(Scenario==2);
CellLine = CellLine_tmpscen(Scenario==2); 
Condition = Condition_tmpscen(Scenario==2);
Day = Day_tmpscen(Scenario==2);
Replicate = Replicate_tmpscen(Scenario==2);
Radius = Radius_tmpscen(Scenario==2);
Outlier = Outlier_tmpscen(Scenario==2);

save([filepath_save 'WM793b_confocal_scenario_2_necrotic' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


fname = fname_tmpscen(Scenario==3);
ExpID = ExpID_tmpscen(Scenario==3);
CellLine = CellLine_tmpscen(Scenario==3); 
Condition = Condition_tmpscen(Scenario==3);
Day = Day_tmpscen(Scenario==3);
Replicate = Replicate_tmpscen(Scenario==3);
Radius = Radius_tmpscen(Scenario==3);
Outlier = Outlier_tmpscen(Scenario==3);

save([filepath_save 'WM793b_confocal_scenario_3_necrotic' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


fname = fname_tmpscen(Scenario==4);
ExpID = ExpID_tmpscen(Scenario==4);
CellLine = CellLine_tmpscen(Scenario==4); 
Condition = Condition_tmpscen(Scenario==4);
Day = Day_tmpscen(Scenario==4);
Replicate = Replicate_tmpscen(Scenario==4);
Radius = Radius_tmpscen(Scenario==4);
Outlier = Outlier_tmpscen(Scenario==4);

save([filepath_save 'WM793b_confocal_scenario_4_necrotic' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


fname = fname_tmpscen(Scenario==5);
ExpID = ExpID_tmpscen(Scenario==5);
CellLine = CellLine_tmpscen(Scenario==5); 
Condition = Condition_tmpscen(Scenario==5);
Day = Day_tmpscen(Scenario==5);
Replicate = Replicate_tmpscen(Scenario==5);
Radius = Radius_tmpscen(Scenario==5);
Outlier = Outlier_tmpscen(Scenario==5);

save([filepath_save 'WM793b_confocal_scenario_5_necrotic' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


%% 7.4) Save PIM outer radius data

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.PIMOuterRadius;
Outlier = zeros(length(confocal_WM793b.AllData.PIMOuterRadius),1);

expected_number_of_outliers = length(Radius) - length(scenario_1_plotdata_PIMOuter) - length(scenario_2_plotdata_PIMOuter);

loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_1_conditions)>0
        % if in list keep
        Scenario(i) = 1;
        if sum((Day(i)== scenario_1_days_to_plot_PIMOuter).*(Radius(i)==scenario_1_plotdata_PIMOuter)) > 0
            Outlier(i) = 0;
        else
            Outlier(i) = 1;
        end
    elseif  sum(Condition(i) == scenario_2_conditions)>0
        Scenario(i) = 2;
        if sum((Day(i)== scenario_2_days_to_plot_PIMOuter).*(Radius(i)==scenario_2_plotdata_PIMOuter)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    elseif  sum(Condition(i) == scenario_3_conditions)>0
        Scenario(i) = 3;
        if sum((Day(i)== scenario_3_days_to_plot_PIMOuter).*(Radius(i)==scenario_3_plotdata_PIMOuter)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    elseif  sum(Condition(i) == scenario_4_conditions)>0
        Scenario(i) = 4;
        if sum((Day(i)== scenario_4_days_to_plot_PIMOuter).*(Radius(i)==scenario_4_plotdata_PIMOuter)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    elseif  sum(Condition(i) == scenario_5_conditions)>0
        Scenario(i) = 5;
        if sum((Day(i)== scenario_5_days_to_plot_PIMOuter).*(Radius(i)==scenario_5_plotdata_PIMOuter)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end

% sum(Outlier)


save([filepath_save 'WM793b_confocal_PIMOuter' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier');


fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;



fname = fname_tmpscen(Scenario==1);
ExpID = ExpID_tmpscen(Scenario==1);
CellLine = CellLine_tmpscen(Scenario==1); 
Condition = Condition_tmpscen(Scenario==1);
Day = Day_tmpscen(Scenario==1);
Replicate = Replicate_tmpscen(Scenario==1);
Radius = Radius_tmpscen(Scenario==1);
Outlier = Outlier_tmpscen(Scenario==1);

save([filepath_save 'WM793b_confocal_scenario_1_PIMOuter' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');

fname = fname_tmpscen(Scenario==2);
ExpID = ExpID_tmpscen(Scenario==2);
CellLine = CellLine_tmpscen(Scenario==2); 
Condition = Condition_tmpscen(Scenario==2);
Day = Day_tmpscen(Scenario==2);
Replicate = Replicate_tmpscen(Scenario==2);
Radius = Radius_tmpscen(Scenario==2);
Outlier = Outlier_tmpscen(Scenario==2);

save([filepath_save 'WM793b_confocal_scenario_2_PIMOuter' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');




fname = fname_tmpscen(Scenario==3);
ExpID = ExpID_tmpscen(Scenario==3);
CellLine = CellLine_tmpscen(Scenario==3); 
Condition = Condition_tmpscen(Scenario==3);
Day = Day_tmpscen(Scenario==3);
Replicate = Replicate_tmpscen(Scenario==3);
Radius = Radius_tmpscen(Scenario==3);
Outlier = Outlier_tmpscen(Scenario==3);

save([filepath_save 'WM793b_confocal_scenario_3_PIMOuter' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


fname = fname_tmpscen(Scenario==4);
ExpID = ExpID_tmpscen(Scenario==4);
CellLine = CellLine_tmpscen(Scenario==4); 
Condition = Condition_tmpscen(Scenario==4);
Day = Day_tmpscen(Scenario==4);
Replicate = Replicate_tmpscen(Scenario==4);
Radius = Radius_tmpscen(Scenario==4);
Outlier = Outlier_tmpscen(Scenario==4);

save([filepath_save 'WM793b_confocal_scenario_4_PIMOuter' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


fname = fname_tmpscen(Scenario==5);
ExpID = ExpID_tmpscen(Scenario==5);
CellLine = CellLine_tmpscen(Scenario==5); 
Condition = Condition_tmpscen(Scenario==5);
Day = Day_tmpscen(Scenario==5);
Replicate = Replicate_tmpscen(Scenario==5);
Radius = Radius_tmpscen(Scenario==5);
Outlier = Outlier_tmpscen(Scenario==5);

save([filepath_save 'WM793b_confocal_scenario_5_PIMOuter' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


%% 7.5) Save data for Scenario 3 - for MCMC

%% 7.5.1) Outer radius

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.OuterRadius;
Outlier = zeros(length(confocal_WM793b.AllData.OuterRadius),1);


loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_3_conditions)>0
        Scenario(i) = 3;
        if sum((Day(i)== scenario_3_days_to_plot_outer).*(Radius(i)==scenario_3_plotdata_outer)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end

fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;


fname = fname_tmpscen(Scenario==3);
ExpID = ExpID_tmpscen(Scenario==3);
CellLine = CellLine_tmpscen(Scenario==3);
Condition = Condition_tmpscen(Scenario==3);
Day = Day_tmpscen(Scenario==3);
Replicate = Replicate_tmpscen(Scenario==3);
Radius = Radius_tmpscen(Scenario==3);
Outlier = Outlier_tmpscen(Scenario==3);

save([filepath_save 'WM793b_confocal_scenario_3_outer' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


% 7.5.2) Save inhibited radius data

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.InhibitedRadius;
Outlier = zeros(length(confocal_WM793b.AllData.InhibitedRadius),1);

loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_3_conditions)>0
        if sum((Day(i)== scenario_3_days_to_plot_inhibited).*(Radius(i)==scenario_3_plotdata_inhibited)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end

fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;


fname = fname_tmpscen(Scenario==3);
ExpID = ExpID_tmpscen(Scenario==3);
CellLine = CellLine_tmpscen(Scenario==3);
Condition = Condition_tmpscen(Scenario==3);
Day = Day_tmpscen(Scenario==3);
Replicate = Replicate_tmpscen(Scenario==3);
Radius = Radius_tmpscen(Scenario==3);
Outlier = Outlier_tmpscen(Scenario==3);

save([filepath_save 'WM793b_confocal_scenario_3_inhibited' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');

% 7.5.3) Save necrotic radius data

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.NecroticRadius;
Outlier = zeros(length(confocal_WM793b.AllData.NecroticRadius),1);


loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_3_conditions)>0
        if sum((Day(i)== scenario_3_days_to_plot_necrotic).*(Radius(i)==scenario_3_plotdata_necrotic)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end


fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;

fname = fname_tmpscen(Scenario==3);
ExpID = ExpID_tmpscen(Scenario==3);
CellLine = CellLine_tmpscen(Scenario==3);
Condition = Condition_tmpscen(Scenario==3);
Day = Day_tmpscen(Scenario==3);
Replicate = Replicate_tmpscen(Scenario==3);
Radius = Radius_tmpscen(Scenario==3);
Outlier = Outlier_tmpscen(Scenario==3);

save([filepath_save 'WM793b_confocal_scenario_3_necrotic' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');



% 7.5.4) Save PIM outer radius data

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.PIMOuterRadius;
Outlier = zeros(length(confocal_WM793b.AllData.PIMOuterRadius),1);


loop_include_counter=0;
for i=1:length(Radius)
    if  sum(Condition(i) == scenario_3_conditions)>0
        if sum((Day(i)== scenario_3_days_to_plot_PIMOuter).*(Radius(i)==scenario_3_plotdata_PIMOuter)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end

fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;

fname = fname_tmpscen(Scenario==3);
ExpID = ExpID_tmpscen(Scenario==3);
CellLine = CellLine_tmpscen(Scenario==3);
Condition = Condition_tmpscen(Scenario==3);
Day = Day_tmpscen(Scenario==3);
Replicate = Replicate_tmpscen(Scenario==3);
Radius = Radius_tmpscen(Scenario==3);
Outlier = Outlier_tmpscen(Scenario==3);

save([filepath_save 'WM793b_confocal_scenario_3_PIMOuter' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


%% 7.6) Save data for Scenario 4 - for MCMC

%% 7.6.1) Outer radius

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.OuterRadius;
Outlier = zeros(length(confocal_WM793b.AllData.OuterRadius),1);

loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_4_conditions)>0
        Scenario(i) = 4;
        if sum((Day(i)== scenario_4_days_to_plot_outer).*(Radius(i)==scenario_4_plotdata_outer)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end

fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;


fname = fname_tmpscen(Scenario==4);
ExpID = ExpID_tmpscen(Scenario==4);
CellLine = CellLine_tmpscen(Scenario==4);
Condition = Condition_tmpscen(Scenario==4);
Day = Day_tmpscen(Scenario==4);
Replicate = Replicate_tmpscen(Scenario==4);
Radius = Radius_tmpscen(Scenario==4);
Outlier = Outlier_tmpscen(Scenario==4);

save([filepath_save 'WM793b_confocal_scenario_4_outer' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


% 7.6.2) Save inhibited radius data

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.InhibitedRadius;
Outlier = zeros(length(confocal_WM793b.AllData.InhibitedRadius),1);


loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_4_conditions)>0
        if sum((Day(i)== scenario_4_days_to_plot_inhibited).*(Radius(i)==scenario_4_plotdata_inhibited)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end

fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;


fname = fname_tmpscen(Scenario==4);
ExpID = ExpID_tmpscen(Scenario==4);
CellLine = CellLine_tmpscen(Scenario==4);
Condition = Condition_tmpscen(Scenario==4);
Day = Day_tmpscen(Scenario==4);
Replicate = Replicate_tmpscen(Scenario==4);
Radius = Radius_tmpscen(Scenario==4);
Outlier = Outlier_tmpscen(Scenario==4);

save([filepath_save 'WM793b_confocal_scenario_4_inhibited' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');

% 7.6.3) Save necrotic radius data

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.NecroticRadius;
Outlier = zeros(length(confocal_WM793b.AllData.NecroticRadius),1);


loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_4_conditions)>0
        if sum((Day(i)== scenario_4_days_to_plot_necrotic).*(Radius(i)==scenario_4_plotdata_necrotic)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end


fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;

fname = fname_tmpscen(Scenario==4);
ExpID = ExpID_tmpscen(Scenario==4);
CellLine = CellLine_tmpscen(Scenario==4);
Condition = Condition_tmpscen(Scenario==4);
Day = Day_tmpscen(Scenario==4);
Replicate = Replicate_tmpscen(Scenario==4);
Radius = Radius_tmpscen(Scenario==4);
Outlier = Outlier_tmpscen(Scenario==4);

save([filepath_save 'WM793b_confocal_scenario_4_necrotic' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');



% 7.6.4) Save PIM outer radius data

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.PIMOuterRadius;
Outlier = zeros(length(confocal_WM793b.AllData.PIMOuterRadius),1);


loop_include_counter=0;
for i=1:length(Radius)
    if  sum(Condition(i) == scenario_4_conditions)>0
        if sum((Day(i)== scenario_4_days_to_plot_PIMOuter).*(Radius(i)==scenario_4_plotdata_PIMOuter)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end

fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;

fname = fname_tmpscen(Scenario==4);
ExpID = ExpID_tmpscen(Scenario==4);
CellLine = CellLine_tmpscen(Scenario==4);
Condition = Condition_tmpscen(Scenario==4);
Day = Day_tmpscen(Scenario==4);
Replicate = Replicate_tmpscen(Scenario==4);
Radius = Radius_tmpscen(Scenario==4);
Outlier = Outlier_tmpscen(Scenario==4);

save([filepath_save 'WM793b_confocal_scenario_4_PIMOuter' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');
	
%% 7.7) Save data for Scenario 5 - for MCMC

%% 7.7.1) Outer radius

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.OuterRadius;
Outlier = zeros(length(confocal_WM793b.AllData.OuterRadius),1);


loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_5_conditions)>0
        Scenario(i) = 5;
        if sum((Day(i)== scenario_5_days_to_plot_outer).*(Radius(i)==scenario_5_plotdata_outer)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end

fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;


fname = fname_tmpscen(Scenario==5);
ExpID = ExpID_tmpscen(Scenario==5);
CellLine = CellLine_tmpscen(Scenario==5);
Condition = Condition_tmpscen(Scenario==5);
Day = Day_tmpscen(Scenario==5);
Replicate = Replicate_tmpscen(Scenario==5);
Radius = Radius_tmpscen(Scenario==5);
Outlier = Outlier_tmpscen(Scenario==5);

save([filepath_save 'WM793b_confocal_scenario_5_outer' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


% 7.7.2) Save inhibited radius data

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.InhibitedRadius;
Outlier = zeros(length(confocal_WM793b.AllData.InhibitedRadius),1);


loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_5_conditions)>0
        if sum((Day(i)== scenario_5_days_to_plot_inhibited).*(Radius(i)==scenario_5_plotdata_inhibited)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end

fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;


fname = fname_tmpscen(Scenario==5);
ExpID = ExpID_tmpscen(Scenario==5);
CellLine = CellLine_tmpscen(Scenario==5);
Condition = Condition_tmpscen(Scenario==5);
Day = Day_tmpscen(Scenario==5);
Replicate = Replicate_tmpscen(Scenario==5);
Radius = Radius_tmpscen(Scenario==5);
Outlier = Outlier_tmpscen(Scenario==5);

save([filepath_save 'WM793b_confocal_scenario_5_inhibited' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');

% 7.7.3) Save necrotic radius data

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.NecroticRadius;
Outlier = zeros(length(confocal_WM793b.AllData.NecroticRadius),1);

loop_include_counter=0;
for i=1:length(Radius)
    if sum(Condition(i) == scenario_5_conditions)>0
        if sum((Day(i)== scenario_5_days_to_plot_necrotic).*(Radius(i)==scenario_5_plotdata_necrotic)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end


fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;

fname = fname_tmpscen(Scenario==5);
ExpID = ExpID_tmpscen(Scenario==5);
CellLine = CellLine_tmpscen(Scenario==5);
Condition = Condition_tmpscen(Scenario==5);
Day = Day_tmpscen(Scenario==5);
Replicate = Replicate_tmpscen(Scenario==5);
Radius = Radius_tmpscen(Scenario==5);
Outlier = Outlier_tmpscen(Scenario==5);

save([filepath_save 'WM793b_confocal_scenario_5_necrotic' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');



% 7.7.4) Save PIM outer radius data

fname = confocal_WM793b.AllData.fname;
ExpID = confocal_WM793b.AllData.ExpID;
CellLine = confocal_WM793b.AllData.CellLine;
Condition = confocal_WM793b.AllData.Condition;
Day = confocal_WM793b.AllData.Day;
Replicate = confocal_WM793b.AllData.Replicate;
Radius = confocal_WM793b.AllData.PIMOuterRadius;
Outlier = zeros(length(confocal_WM793b.AllData.PIMOuterRadius),1);


loop_include_counter=0;
for i=1:length(Radius)
    if  sum(Condition(i) == scenario_5_conditions)>0
        if sum((Day(i)== scenario_5_days_to_plot_PIMOuter).*(Radius(i)==scenario_5_plotdata_PIMOuter)) > 0
            Outlier(i) = 0;
        else
            Outlier(i)=1;
        end
    end
end

fname_tmpscen = fname;
ExpID_tmpscen = ExpID;
CellLine_tmpscen = CellLine;
Condition_tmpscen = Condition;
Day_tmpscen = Day;
Replicate_tmpscen = Replicate;
Radius_tmpscen = Radius;
Outlier_tmpscen =  Outlier;

fname = fname_tmpscen(Scenario==5);
ExpID = ExpID_tmpscen(Scenario==5);
CellLine = CellLine_tmpscen(Scenario==5);
Condition = Condition_tmpscen(Scenario==5);
Day = Day_tmpscen(Scenario==5);
Replicate = Replicate_tmpscen(Scenario==5);
Radius = Radius_tmpscen(Scenario==5);
Outlier = Outlier_tmpscen(Scenario==5);

save([filepath_save 'WM793b_confocal_scenario_5_PIMOuter' '.mat'],'-v7.3',...
    'fname',...
    'ExpID',...
    'CellLine',...
    'Condition',...
    'Day',...
    'Replicate',...
    'Radius',...
    'Outlier',...
    'Scenario');


%% 8) Change directory back to 2Code
cd ../

