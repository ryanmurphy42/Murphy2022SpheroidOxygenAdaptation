%% WM164 oxygen and waste analysis
% Using mathematical modelling and statistical analysis to interpret the experimental data
% Normoxia vs Hypoxia
% Spheroid snapshots
% Estimate Rc, alpha, p_i, Rd (mathcalR)

% 1.1) Load data
% 2) Estimate R_c, alpha and validate with PIMOuter measurements
% 3) PLOTS - Rcrit
% 4) PLOTS - alpha
% 5) PLOTS - p_n (should be zero)
% 6) PLOTS - Measured Rp(t) v Predicted Rp(t)
% 7) PLOT - Rd
% 8) PLOT - p_i
% 9) Change directory back to 2Code


%% 1.1) Load data
% open data from 1Data folder
mydir  = pwd; idcs   = strfind(mydir,'\'); newdir = mydir(1:idcs(end)-1);

confocal.AllData = readtable([newdir '\1Data\ExpOA_WM164.xlsx'],'Sheet', 'Sheet1'); % load data
wm164_data_include = (sum(cell2mat(confocal.AllData.CellLine) == 'WM164',2)==5); % identify WM164 data
confocal_WM164.AllDatatmp2 = confocal.AllData(wm164_data_include,:); % select WM164 data
confocal_WM164.AllDatatmp1 = confocal_WM164.AllDatatmp2((confocal_WM164.AllDatatmp2.Keep==1),:); % only keep those that are not multiple spheroids or damaged spheroids
confocal_WM164.AllData = confocal_WM164.AllDatatmp1(sum(confocal_WM164.AllDatatmp1.Condition==[1,2,3,4,17,18,19,15,16],2)==1,:); % only include scenarios 1 and 2

% create folders
cd ../
if ~exist('3MATLABsavefiles\1FirstAnalysisExperimentalData\WM164', 'dir')
    mkdir(['3MATLABsavefiles\1FirstAnalysisExperimentalData\WM164']);
end
if ~exist('4Figures\1FirstAnalysisExperimentalData\WM164', 'dir')
    mkdir(['4Figures\1FirstAnalysisExperimentalData\WM164']);
end
if ~exist('4Figures\1FirstAnalysisExperimentalData\WM164\OxygenProfiles', 'dir')
    mkdir(['4Figures\1FirstAnalysisExperimentalData\WM164\OxygenProfiles']);
end

cd '2Code';
cd '1FirstAnalysisExperimentalData';

filepath_save = [newdir '\3MATLABsavefiles\1FirstAnalysisExperimentalData\WM164\'];
filepath_save_figs =  [newdir '\4Figures\1FirstAnalysisExperimentalData\WM164\'];

% setting for figures plots
jitter_plot = 0.18;
map_custom = [0 0 1
    1 0 0];
include_oxygen_profile_plots = 0;


%% 2) Estimate R_c, alpha and validate with PIMOuter measurements

Omega_unconverted = 3.0318*10^7; %mmHg kg m^(-3)
Omega = Omega_unconverted*(21/160)*(1*10^(-6))^3; % [% kg micrometers^(-3)]

p_infinity_21 = 21; % [percent]
p_infinity_2 = 2; % [percent]

k_unconverted = 2*10^(-9); % [m^2 s^-1]
k = k_unconverted*(1/(1*10^(-6))^2)*(1./(60*60*24)); % [micrometre^2 day^(-1)]

% Compute Rc, alpha, pi, Rd for each spheroid

number_of_spheroids = size(confocal_WM164.AllData,1);

Rcrit_vector = zeros(number_of_spheroids,1); %[micrometres]
alpha_vector = zeros(number_of_spheroids,1); %[m^3 kg^-1 s-1]
p_n_vector= zeros(number_of_spheroids,1); %[%]
p_i_vector = zeros(number_of_spheroids,1); %[%]
PIMOuterRadius_prediction_vector  = zeros(number_of_spheroids,1); %[micrometres]
Rd_vector = zeros(number_of_spheroids,1); %[micrometres]

p_pim_activate = 1.32;

scenario_colour_vec = zeros(number_of_spheroids,1);

for i=1:number_of_spheroids

    if  sum(confocal_WM164.AllData.Condition(i) == [5,6,7,8,9,10,15,16,17,18,19])> 0
        % Measurement at 21% oxygen
        p_infinity = p_infinity_21;
    elseif sum(confocal_WM164.AllData.Condition(i) == [1,2,3,4,11,12,13,14]) > 0
        % Measurement at 2% oxygen
        p_infinity = p_infinity_2;
    end

    if  sum(confocal_WM164.AllData.Condition(i) == [1,2,3,4])> 0
        % Scenario 1
        scenario_colour_vec(i) = 0;
    elseif sum(confocal_WM164.AllData.Condition(i) == [17,18,19,15,16]) > 0
        % Scenario 2
        scenario_colour_vec(i) = 1;
    else
        % other scenario
        scenario_colour_vec(i) = -1;
    end

    OuterRadius_this_loop = confocal_WM164.AllData.OuterRadius(i);
    NecroticRadius_this_loop = confocal_WM164.AllData.NecroticRadius(i);
    PIMOuterRadius_this_loop = confocal_WM164.AllData.PIMOuterRadius(i);
    InhibitedRadius_this_loop = confocal_WM164.AllData.InhibitedRadius(i);


    % for each spheroid compute R_c [micrometres]
    Rcrit_this_loop = sqrt((OuterRadius_this_loop.^2 - NecroticRadius_this_loop.^2) - (2*NecroticRadius_this_loop.^2/OuterRadius_this_loop)*(OuterRadius_this_loop - NecroticRadius_this_loop));

    % for each spheroid compute alpha [micrometre^3 kg^-1 day^-1]
    alpha_this_loop = (6*k/(Omega*Rcrit_this_loop.^2))*p_infinity;

    % for each spheroid compute alpha [m^3 kg^-1 day^-1]
    alpha_this_loop_m3kgminus1secondsminus1 = alpha_this_loop/(1/(1*10^(-6))^3)/(1./(60*60*24));


    % for each spheroid compute p_n [%]
    p_n_this_loop =   p_infinity - (alpha_this_loop.*Omega./(6*k)).*(OuterRadius_this_loop.^2 - NecroticRadius_this_loop.^2) + ...
        (alpha_this_loop.*Omega./(3*k)).*((NecroticRadius_this_loop.^2)./OuterRadius_this_loop).*(OuterRadius_this_loop - NecroticRadius_this_loop);


    % for each spheroid compute Rp [micrometres] = radial position where equal to p_pim =1.32
    function_oxygen_partial_pressure = @(r) + p_infinity - (alpha_this_loop.*Omega./(6*k)).*(OuterRadius_this_loop.^2 - r.^2) + ...
        (alpha_this_loop.*Omega./(3*k)).*(NecroticRadius_this_loop.^3).*(1./r - 1./OuterRadius_this_loop);

    PIMOuterRadius_prediction_this_loop = fzero(@(r) function_oxygen_partial_pressure(r) - p_pim_activate ...
        ,[NecroticRadius_this_loop+.01,OuterRadius_this_loop]);



    % estimate p_i [%]
    p_i_this_loop = -1;
    if InhibitedRadius_this_loop > 0 % if in phase 2 or 3
        if InhibitedRadius_this_loop > NecroticRadius_this_loop %
            p_i_this_loop = p_infinity - (alpha_this_loop.*Omega./(6*k)).*(OuterRadius_this_loop.^2 - InhibitedRadius_this_loop.^2) + ...
                (alpha_this_loop.*Omega./(3*k)).*(NecroticRadius_this_loop.^3).*(1./InhibitedRadius_this_loop - 1./OuterRadius_this_loop);
        end
    end

    % oxygen profile figure
    if include_oxygen_profile_plots == 1
        r_plot_1 = linspace(NecroticRadius_this_loop,OuterRadius_this_loop,100);
        r_plot_2 = linspace(0,NecroticRadius_this_loop,2);
        figure
        hold on
        plot([PIMOuterRadius_this_loop,PIMOuterRadius_this_loop],[0,p_pim_activate],'LineWidth',3,'Color',	'#00FFFF') % actual PIM
        plot([PIMOuterRadius_prediction_this_loop,PIMOuterRadius_prediction_this_loop],[0,p_pim_activate],'LineWidth',3,'Color','#7E2F8E') % predicted PIM
        plot([InhibitedRadius_this_loop,InhibitedRadius_this_loop],[0,p_i_this_loop],'LineWidth',3,'Color',		'm') % inhibited radius
        plot(r_plot_1,p_infinity - (alpha_this_loop.*Omega./(6*k)).*(OuterRadius_this_loop.^2 - r_plot_1.^2) + ...
            (alpha_this_loop.*Omega./(3*k)).*(NecroticRadius_this_loop.^3).*(1./r_plot_1 - 1./OuterRadius_this_loop),'LineWidth',3,'Color','#EDB120'); % oxygen profile
        plot(r_plot_2,[0,0],'LineWidth',3,'Color',		'#D95319')
        legend('Actual PIM', 'Prediction PIM', 'Inhibited')
        % necrotic core boundary
        xlabel('Radius')
        ylabel('Oxygen [%]')
        grid on
        box on
        ylim([0,25])
        xlim([0,500])
        saveas(gcf, [filepath_save_figs '\OxygenProfiles\' 'OxygenProfile_' num2str(i) '_' char(confocal_WM164.AllData.fname(i)) '.fig'])
        saveas(gcf, [filepath_save_figs '\OxygenProfiles\' 'OxygenProfile_' num2str(i) '_' char(confocal_WM164.AllData.fname(i)) '.png'])
        print([filepath_save_figs '\OxygenProfiles\' 'OxygenProfile_'  num2str(i) '_' char(confocal_WM164.AllData.fname(i))],'-depsc2','-painters')
        close
    end



    Rdsquared_this_loop = -1;
    Rd_this_loop = -1;
    if InhibitedRadius_this_loop > 0
        if InhibitedRadius_this_loop > NecroticRadius_this_loop
            % estimate Rdsquared [micrometres^2]
            Rdsquared_this_loop = (OuterRadius_this_loop.^2 - InhibitedRadius_this_loop.^2) - 2.*(NecroticRadius_this_loop.^3).*(1./InhibitedRadius_this_loop - 1./OuterRadius_this_loop);
            % estimate Rd [micrometres]
            Rd_this_loop = sqrt(Rdsquared_this_loop);
        end
    end


    % save to vector
    Rcrit_vector(i) = Rcrit_this_loop; % [micrometres]
    alpha_vector(i) = alpha_this_loop_m3kgminus1secondsminus1; % [m^3 kg^-1 s^-1]
    p_n_vector(i) = p_n_this_loop; %[%]
    p_i_vector(i) = p_i_this_loop; %[%]
    PIMOuterRadius_prediction_vector(i) = PIMOuterRadius_prediction_this_loop;% [micrometres]
    Rd_vector(i) = Rd_this_loop; % [micrometres]

end



%% 3) PLOTS - Rcrit


figure
scatter(1:number_of_spheroids,Rcrit_vector,25,confocal_WM164.AllData.Condition,'filled')
xlabel('Spheroid')
ylabel('Rcrit')
ylim([0,500])
hold on
plot([0,number_of_spheroids],[mean(Rcrit_vector(scenario_colour_vec==1)),mean(Rcrit_vector(scenario_colour_vec==1))],'r--')
plot([0,number_of_spheroids],[mean(Rcrit_vector(scenario_colour_vec==0)),mean(Rcrit_vector(scenario_colour_vec==0))],'b--')
box on
grid on
colormap winter
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_1'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_1'],'-depsc2','-painters')


figure
scatter(confocal_WM164.AllData.OuterRadius,Rcrit_vector,25,confocal_WM164.AllData.Condition,'filled')
xlabel('Outer Radius')
ylabel('Rcrit')
xlim([0,500])
ylim([0,500])
box on
grid on
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_2'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_2'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_2'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_2'],'-depsc2','-painters')


figure
scatter(confocal_WM164.AllData.Day,Rcrit_vector,25,scenario_colour_vec,'filled')
hold on
plot([0,8.5],[mean(Rcrit_vector(scenario_colour_vec==1)),mean(Rcrit_vector(scenario_colour_vec==1))],'r--')
plot([0,8.5],[mean(Rcrit_vector(scenario_colour_vec==0)),mean(Rcrit_vector(scenario_colour_vec==0))],'b--')
xlabel('time [days]')
ylabel('Rcrit')
ylim([0,500])
xlim([0,8.5])
box on
grid on
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_4'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_4'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_4'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_4'],'-depsc2','-painters')


% Rcrit box plot - comparing normoxia v hypoxia

Rcrit_vector_normoxia = Rcrit_vector(scenario_colour_vec==1);
Rcrit_vector_hypoxia = Rcrit_vector(scenario_colour_vec==0);

figure
hold on
boxchart(scenario_colour_vec(scenario_colour_vec==1)-1,Rcrit_vector_normoxia,'BoxFaceColor','r','JitterOutliers', 'on','MarkerColor','r')
boxchart(scenario_colour_vec(scenario_colour_vec==0)+1,Rcrit_vector_hypoxia,'BoxFaceColor','b','JitterOutliers', 'on','MarkerColor','b')
% exclude outliers from scatter
%use isoutlier function
Rcrit_indicator_vector_normoxia_to_plot_tmp = scenario_colour_vec(scenario_colour_vec==1);
Rcrit_vector_normoxia_to_plot_tmp =Rcrit_vector_normoxia;
Rcrit_indicator_vector_hypoxia_to_plot_tmp =scenario_colour_vec(scenario_colour_vec==0);
Rcrit_vector_hypoxia_to_plot_tmp =Rcrit_vector_hypoxia;
[Rcrit_indicator_vector_normoxia_to_plot,Rcrit_vector_normoxia_to_plot] = function_removeoutliers_for_plot(Rcrit_indicator_vector_normoxia_to_plot_tmp,Rcrit_vector_normoxia_to_plot_tmp);
[Rcrit_indicator_vector_hypoxia_to_plot,Rcrit_vector_hypoxia_to_plot] = function_removeoutliers_for_plot(Rcrit_indicator_vector_hypoxia_to_plot_tmp,Rcrit_vector_hypoxia_to_plot_tmp);
s2 = scatter(Rcrit_indicator_vector_normoxia_to_plot-1,Rcrit_vector_normoxia_to_plot,'r','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(Rcrit_indicator_vector_hypoxia_to_plot+1,Rcrit_vector_hypoxia_to_plot,'b','filled','jitter', 'on', 'jitterAmount', jitter_plot);
alpha(s2,0.2)
alpha(s1,0.2)
xlim([-0.5,1.5])
ylim([0,500])
xlabel('Normoxia | Hypoxia')
ylabel('Rcrit [microns]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_10'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_10'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_10'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rcrit_10'],'-depsc2','-painters')

%% 4) PLOTS - alpha


figure
scatter(1:number_of_spheroids,alpha_vector,25,confocal_WM164.AllData.Condition,'filled')
xlabel('Spheroid')
ylabel('alpha')
ylim([0,2.5*10^-6])
box on
grid on
colormap winter
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_1'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_1'],'-depsc2','-painters')


figure
scatter(confocal_WM164.AllData.OuterRadius,alpha_vector,25,confocal_WM164.AllData.Condition,'filled')
xlabel('Outer Radius')
ylabel('alpha')
xlim([0,500])
ylim([0,2.5*10^-6])
box on
grid on
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_2'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_2'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_2'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_2'],'-depsc2','-painters')


figure
scatter(confocal_WM164.AllData.Day,alpha_vector,25,scenario_colour_vec,'filled')
xlabel('time [days]')
ylabel('alpha')
ylim([0,2.5*10^-6])
xlim([0,8.5])
box on
grid on
colormap jet
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_4'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_4'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_4'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_4'],'-depsc2','-painters')



alpha_vector_normoxia = alpha_vector(scenario_colour_vec==1);
alpha_vector_hypoxia = alpha_vector(scenario_colour_vec==0);

figure
hold on
boxchart(scenario_colour_vec(scenario_colour_vec==1)-1,alpha_vector_normoxia,'BoxFaceColor','r','JitterOutliers', 'on','MarkerColor','r')
boxchart(scenario_colour_vec(scenario_colour_vec==0)+1,alpha_vector_hypoxia,'BoxFaceColor','b','JitterOutliers', 'on','MarkerColor','b')
% exclude outliers from scatter
%use isoutlier function
alpha_indicator_vector_normoxia_to_plot_tmp = scenario_colour_vec(scenario_colour_vec==1);
alpha_vector_normoxia_to_plot_tmp =alpha_vector_normoxia;
alpha_indicator_vector_hypoxia_to_plot_tmp =scenario_colour_vec(scenario_colour_vec==0);
alpha_vector_hypoxia_to_plot_tmp =alpha_vector_hypoxia;
[alpha_indicator_vector_normoxia_to_plot,alpha_vector_normoxia_to_plot] = function_removeoutliers_for_plot(alpha_indicator_vector_normoxia_to_plot_tmp,alpha_vector_normoxia_to_plot_tmp);
[alpha_indicator_vector_hypoxia_to_plot,alpha_vector_hypoxia_to_plot] = function_removeoutliers_for_plot(alpha_indicator_vector_hypoxia_to_plot_tmp,alpha_vector_hypoxia_to_plot_tmp);
s2 = scatter(alpha_indicator_vector_normoxia_to_plot-1,alpha_vector_normoxia_to_plot,'r','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(alpha_indicator_vector_hypoxia_to_plot+1,alpha_vector_hypoxia_to_plot,'b','filled','jitter', 'on', 'jitterAmount', jitter_plot);
alpha(s2,0.2)
alpha(s1,0.2)
xlim([-0.5,1.5])
ylim([0,2.5*10^-6])
xlabel('Normoxia | Hypoxia')
ylabel('alpha [units]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_10'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_10'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_10'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_alpha_10'],'-depsc2','-painters')


%% 5) PLOTS - p_n (should be zero)



figure
scatter(1:number_of_spheroids,p_n_vector,25,confocal_WM164.AllData.Condition,'filled')
xlabel('Spheroid')
ylabel('pn')
% ylim([0,500])
box on
grid on
colormap winter
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pn_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pn_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pn_1'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pn_1'],'-depsc2','-painters')




%% 6) PLOTS - Measured Rp(t) v Predicted Rp(t)




figure
hold on
scatter(confocal_WM164.AllData.PIMOuterRadius(confocal_WM164.AllData.PIMOuterRadius>0),PIMOuterRadius_prediction_vector(confocal_WM164.AllData.PIMOuterRadius>0),25,scenario_colour_vec(confocal_WM164.AllData.PIMOuterRadius>0),'filled')
xlabel('Measured OuterPIM')
ylabel('Predicted OuterPIM')
plot([0,500],[0,500],'K--')
xlim([0,500])
ylim([0,500])
box on
grid on
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pim_6'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pim_6'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pim_6'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pim_6'],'-depsc2','-painters')


mdl_PIMOuterRadius = fitlm(confocal_WM164.AllData.PIMOuterRadius(confocal_WM164.AllData.PIMOuterRadius>0),PIMOuterRadius_prediction_vector(confocal_WM164.AllData.PIMOuterRadius>0));
figure
hold on
scatter(confocal_WM164.AllData.PIMOuterRadius(confocal_WM164.AllData.PIMOuterRadius>0),PIMOuterRadius_prediction_vector(confocal_WM164.AllData.PIMOuterRadius>0),25,scenario_colour_vec(confocal_WM164.AllData.PIMOuterRadius>0),'filled')
x_plot_OxygenAnalysis_Scenario_1_2_pim_7=0:500;
y_plot_OxygenAnalysis_Scenario_1_2_pim_7=mdl_PIMOuterRadius.Coefficients.Estimate(1) + mdl_PIMOuterRadius.Coefficients.Estimate(2).*x_plot_OxygenAnalysis_Scenario_1_2_pim_7;
plot(x_plot_OxygenAnalysis_Scenario_1_2_pim_7,y_plot_OxygenAnalysis_Scenario_1_2_pim_7)
xlabel('Measured OuterPIM')
ylabel('Predicted OuterPIM')
plot([0,500],[0,500],'K--')
xlim([0,500])
ylim([0,500])
box on
grid on
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pim_7'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pim_7'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pim_7'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pim_7'],'-depsc2','-painters')


%% 7) PLOT - p_i




figure
hold on
scatter(1:number_of_spheroids,p_i_vector,25,confocal_WM164.AllData.Condition,'filled')
plot([0,number_of_spheroids],[21,21],'k--')
plot([0,number_of_spheroids],[2,2],'k--')
xlabel('Spheroid')
ylabel('pi')
ylim([0,12])
box on
grid on
colormap winter
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_1'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_1'],'-depsc2','-painters')


figure
hold on
scatter(confocal_WM164.AllData.OuterRadius,p_i_vector,25,confocal_WM164.AllData.Condition,'filled')
plot([0,500],[21,21],'k--')
plot([0,500],[2,2],'k--')
xlabel('Outer Radius')
ylabel('pi')
xlim([0,500])
ylim([0,12])
box on
grid on
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_2'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_2'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_2'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_2'],'-depsc2','-painters')


figure
hold on
scatter(confocal_WM164.AllData.Day,p_i_vector,25,scenario_colour_vec,'filled')
plot([0,8.5],[21,21],'k--')
plot([0,8.5],[2,2],'k--')
xlabel('time [days]')
ylabel('pi')
ylim([0,12])
xlim([0,8.5])
box on
grid on
colormap jet
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_4'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_4'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_4'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_4'],'-depsc2','-painters')


figure
hold on
scatter(confocal_WM164.AllData.NecroticRadius,p_i_vector,25,scenario_colour_vec,'filled')
plot([0,500],[21,21],'k--')
plot([0,500],[2,2],'k--')
xlabel('Necrotic Radius')
ylabel('pi')
xlim([0,500])
ylim([0,10])
box on
grid on
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_5'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_5'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_5'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_5'],'-depsc2','-painters')


figure
hold on
scatter(confocal_WM164.AllData.InhibitedRadius,p_i_vector,25,scenario_colour_vec,'filled')
plot([0,500],[21,21],'k--')
plot([0,500],[2,2],'k--')
xlabel('Inhibited Radius')
ylabel('pi')
xlim([0,500])
ylim([0,10])
box on
grid on
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_6'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_6'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_6'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_6'],'-depsc2','-painters')


% pi box plot - comparing normoxia v hypoxia

p_i_vector_tmp_2 = p_i_vector;
p_i_vector_tmp_2 = p_i_vector_tmp_2(p_i_vector_tmp_2>=0);
scenario_colour_vec_tmp_2 = scenario_colour_vec;
scenario_colour_vec_tmp_2 = scenario_colour_vec(p_i_vector_tmp_2>=0);

pi_vector_normoxia = p_i_vector_tmp_2(scenario_colour_vec_tmp_2==1);
pi_vector_hypoxia = p_i_vector_tmp_2(scenario_colour_vec_tmp_2==0);

figure
hold on
plot([-.5,1.5],[21,21],'k--')
plot([-.5,1.5],[2,2],'k--')
boxchart(scenario_colour_vec_tmp_2(scenario_colour_vec_tmp_2==1)-1,pi_vector_normoxia,'BoxFaceColor','r','JitterOutliers', 'on','MarkerColor','r')
boxchart(scenario_colour_vec_tmp_2(scenario_colour_vec_tmp_2==0)+1,pi_vector_hypoxia,'BoxFaceColor','b','JitterOutliers', 'on','MarkerColor','b')
% exclude outliers from scatter
%use isoutlier function
pi_indicator_vector_normoxia_to_plot_tmp = scenario_colour_vec_tmp_2(scenario_colour_vec_tmp_2==1);
pi_vector_normoxia_to_plot_tmp =pi_vector_normoxia;
pi_indicator_vector_hypoxia_to_plot_tmp =scenario_colour_vec_tmp_2(scenario_colour_vec_tmp_2==0);
pi_vector_hypoxia_to_plot_tmp =pi_vector_hypoxia;
[pi_indicator_vector_normoxia_to_plot,pi_vector_normoxia_to_plot] = function_removeoutliers_for_plot(pi_indicator_vector_normoxia_to_plot_tmp,pi_vector_normoxia_to_plot_tmp);
[pi_indicator_vector_hypoxia_to_plot,pi_vector_hypoxia_to_plot] = function_removeoutliers_for_plot(pi_indicator_vector_hypoxia_to_plot_tmp,pi_vector_hypoxia_to_plot_tmp);
s2 = scatter(pi_indicator_vector_normoxia_to_plot-1,pi_vector_normoxia_to_plot,'r','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(pi_indicator_vector_hypoxia_to_plot+1,pi_vector_hypoxia_to_plot,'b','filled','jitter', 'on', 'jitterAmount', jitter_plot);
alpha(s2,0.2)
alpha(s1,0.2)
xlim([-0.5,1.5])
ylim([0,12])
xlabel('Normoxia | Hypoxia')
ylabel('p_i [%]')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_10'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_10'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_10'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_pi_10'],'-depsc2','-painters')



%% 8) PLOT - Rd

figure
hold on
scatter(1:number_of_spheroids,Rd_vector,25,confocal_WM164.AllData.Condition,'filled')
xlabel('Spheroid')
ylabel('Rd')
box on
grid on
colormap(map_custom)
ylim([0,500])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_1'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_1'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_1'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_1'],'-depsc2','-painters')


figure
scatter(confocal_WM164.AllData.OuterRadius,Rd_vector,25,confocal_WM164.AllData.Condition,'filled')
xlabel('Outer Radius')
ylabel('Rd')
xlim([0,500])
ylim([0,500])
box on
grid on
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_2'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_2'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_2'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_2'],'-depsc2','-painters')


figure
scatter(confocal_WM164.AllData.Day,Rd_vector,25,scenario_colour_vec,'filled')
xlabel('time [days]')
ylabel('Rd')
xlim([0,8.5])
ylim([0,500])
box on
grid on
colormap jet
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_4'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_4'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_4'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_4'],'-depsc2','-painters')


figure
hold on
scatter(confocal_WM164.AllData.NecroticRadius,Rd_vector,25,scenario_colour_vec,'filled')
xlabel('Necrotic Radius')
ylabel('Rd')
xlim([0,500])
ylim([0,500])
box on
grid on
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_5'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_5'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_5'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_5'],'-depsc2','-painters')


figure
hold on
scatter(confocal_WM164.AllData.InhibitedRadius,Rd_vector,25,scenario_colour_vec,'filled')
xlabel('Inhibited Radius')
ylabel('Rd')
xlim([0,500])
ylim([0,500])
box on
grid on
colormap(map_custom)
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_6'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_6'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_6'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_6'],'-depsc2','-painters')





% Rd box plot - comparing normoxia v hypoxia

Rd_vector_tmp_2 = Rd_vector;
Rd_vector_tmp_2 = Rd_vector_tmp_2(Rd_vector_tmp_2>=0);
scenario_colour_vec_tmp_2 = scenario_colour_vec;
scenario_colour_vec_tmp_2 = scenario_colour_vec(Rd_vector_tmp_2>=0);

Rd_vector_normoxia = Rd_vector_tmp_2(scenario_colour_vec_tmp_2==1);
Rd_vector_hypoxia = Rd_vector_tmp_2(scenario_colour_vec_tmp_2==0);
figure
hold on
boxchart(scenario_colour_vec_tmp_2(scenario_colour_vec_tmp_2==1)-1,Rd_vector_normoxia,'BoxFaceColor','r','JitterOutliers', 'on','MarkerColor','r')
boxchart(scenario_colour_vec_tmp_2(scenario_colour_vec_tmp_2==0)+1,Rd_vector_hypoxia,'BoxFaceColor','b','JitterOutliers', 'on','MarkerColor','b')
% exclude outliers from scatter
%use isoutlier function
Rd_indicator_vector_normoxia_to_plot_tmp = scenario_colour_vec_tmp_2(scenario_colour_vec_tmp_2==1);
Rd_vector_normoxia_to_plot_tmp =Rd_vector_normoxia;
Rd_indicator_vector_hypoxia_to_plot_tmp =scenario_colour_vec_tmp_2(scenario_colour_vec_tmp_2==0);
Rd_vector_hypoxia_to_plot_tmp =Rd_vector_hypoxia;
[Rd_indicator_vector_normoxia_to_plot,Rd_vector_normoxia_to_plot] = function_removeoutliers_for_plot(Rd_indicator_vector_normoxia_to_plot_tmp,Rd_vector_normoxia_to_plot_tmp);
[Rd_indicator_vector_hypoxia_to_plot,Rd_vector_hypoxia_to_plot] = function_removeoutliers_for_plot(Rd_indicator_vector_hypoxia_to_plot_tmp,Rd_vector_hypoxia_to_plot_tmp);
s2 = scatter(Rd_indicator_vector_normoxia_to_plot-1,Rd_vector_normoxia_to_plot,'r','filled','jitter', 'on', 'jitterAmount', jitter_plot);
s1 = scatter(Rd_indicator_vector_hypoxia_to_plot+1,Rd_vector_hypoxia_to_plot,'b','filled','jitter', 'on', 'jitterAmount', jitter_plot);
alpha(s2,0.2)
alpha(s1,0.2)
xlim([-0.5,1.5])
ylim([0,500])
xlabel('Normoxia | Hypoxia')
ylabel('Rd')
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_10'  '.fig'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_10'  '.pdf'])
saveas(gcf, [filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_10'  '.png'])
print([filepath_save_figs  'OxygenAnalysis_Scenario_1_2_Rd_10'],'-depsc2','-painters')


close all
%% 9) Change directory back to 2Code
cd ../

