function [OuterRadius_to_plot,...
    InhibitedRadius_to_plot,...
    NecroticRadius_to_plot,...
    PIMOuterRadius_to_plot,...
    PIMInnerRadius_to_plot,...
    Days_to_plot,...
    OuterRadius_mean_to_plot,...
    InhibitedRadius_mean_to_plot,...
    NecroticRadius_mean_to_plot,...
    PIMOuterRadius_mean_to_plot,...
    PIMInnerRadius_mean_to_plot,...
    Days_unique_to_plot,...
    OuterRadius_std_to_plot,...
    InhibitedRadius_std_to_plot,...
    NecroticRadius_std_to_plot,...
    PIMOuterRadius_std_to_plot,...
    PIMInnerRadius_std_to_plot,...
    OuterRadius_count_to_plot,...
    InhibitedRadius_count_to_plot,...
    NecroticRadius_count_to_plot,...
    PIMOuterRadius_count_to_plot,...
    PIMInnerRadius_count_to_plot]  =  function_plotdata_scenario(scenario_conditions,...
    Condition,...
    OuterRadius,...
    InhibitedRadius,...
    NecroticRadius,...
    PIMOuterRadius,...
    PIMInnerRadius,...
    Day)

OuterRadius_to_plot = [];
InhibitedRadius_to_plot = [];
NecroticRadius_to_plot = [];
PIMOuterRadius_to_plot = [];
PIMInnerRadius_to_plot = [];
Days_to_plot = [];

OuterRadius_mean_to_plot = [];
InhibitedRadius_mean_to_plot = [];
NecroticRadius_mean_to_plot = [];
PIMOuterRadius_mean_to_plot = [];
PIMInnerRadius_mean_to_plot = [];
Days_unique_to_plot = [];

OuterRadius_std_to_plot = [];
InhibitedRadius_std_to_plot = [];
NecroticRadius_std_to_plot = [];
PIMOuterRadius_std_to_plot = [];
PIMInnerRadius_std_to_plot = [];


OuterRadius_count_to_plot = [];
InhibitedRadius_count_to_plot = [];
NecroticRadius_count_to_plot = [];
PIMOuterRadius_count_to_plot = [];
PIMInnerRadius_count_to_plot = [];


for scenario_condition_counter = 1:length(scenario_conditions)

    % for each scenario select the relevant data
    OuterRadius_tmp = OuterRadius((Condition == scenario_conditions(scenario_condition_counter)),:);
    InhibitedRadius_tmp = InhibitedRadius((Condition == scenario_conditions(scenario_condition_counter)),:);
    NecroticRadius_tmp = NecroticRadius((Condition == scenario_conditions(scenario_condition_counter)),:);
    Days_tmp = Day((Condition == scenario_conditions(scenario_condition_counter)),:);
    PIMOuterRadius_tmp = PIMOuterRadius((Condition == scenario_conditions(scenario_condition_counter)),:);
    PIMInnerRadius_tmp = PIMInnerRadius((Condition == scenario_conditions(scenario_condition_counter)),:);

    % save all data to vectors
    OuterRadius_to_plot = [OuterRadius_to_plot;OuterRadius_tmp];
    InhibitedRadius_to_plot = [InhibitedRadius_to_plot;InhibitedRadius_tmp];
    NecroticRadius_to_plot = [NecroticRadius_to_plot;NecroticRadius_tmp];
    PIMOuterRadius_to_plot = [PIMOuterRadius_to_plot;PIMOuterRadius_tmp];
    PIMInnerRadius_to_plot = [PIMInnerRadius_to_plot;PIMInnerRadius_tmp];
    Days_to_plot = [Days_to_plot;Days_tmp];

    % compute mean and standard deviation of measurements
    OuterRadius_mean_to_plot = [OuterRadius_mean_to_plot;mean(OuterRadius_tmp)];
    InhibitedRadius_mean_to_plot = [InhibitedRadius_mean_to_plot;mean(InhibitedRadius_tmp)];
    NecroticRadius_mean_to_plot = [NecroticRadius_mean_to_plot;mean(NecroticRadius_tmp)];
    PIMOuterRadius_mean_to_plot = [PIMOuterRadius_mean_to_plot;mean(PIMOuterRadius_tmp)];
    PIMInnerRadius_mean_to_plot = [PIMInnerRadius_to_plot;mean(PIMInnerRadius_tmp)];
    Days_unique_to_plot = [Days_unique_to_plot;unique(Days_tmp)];

    OuterRadius_std_to_plot = [OuterRadius_std_to_plot;std(OuterRadius_tmp)];
    InhibitedRadius_std_to_plot = [InhibitedRadius_std_to_plot;std(InhibitedRadius_tmp)];
    NecroticRadius_std_to_plot = [NecroticRadius_std_to_plot;std(NecroticRadius_tmp)];
    PIMOuterRadius_std_to_plot = [PIMOuterRadius_std_to_plot;std(PIMOuterRadius_tmp)];
    PIMInnerRadius_std_to_plot = [PIMInnerRadius_std_to_plot;std(PIMInnerRadius_tmp)];

    % compute total number of measurements
    OuterRadius_count_to_plot = [OuterRadius_count_to_plot;length(OuterRadius_tmp)];
    InhibitedRadius_count_to_plot = [InhibitedRadius_count_to_plot;length(InhibitedRadius_tmp)];
    NecroticRadius_count_to_plot = [NecroticRadius_count_to_plot;length(NecroticRadius_tmp)];
    PIMOuterRadius_count_to_plot = [PIMOuterRadius_count_to_plot;length(PIMOuterRadius_tmp)];
    PIMInnerRadius_count_to_plot = [PIMInnerRadius_count_to_plot;length(PIMInnerRadius_tmp)];

end


end