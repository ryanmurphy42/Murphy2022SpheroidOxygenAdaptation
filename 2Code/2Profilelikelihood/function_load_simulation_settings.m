function  [simulation_id, ...
    data_file_to_use_array,...
    data_sheet_to_use_array,...
    p_lower_bounds,...
    p_upper_bounds,...
    p_first_guess,...
    times_to_use,...
    initial_condition,...
    filepath_save,...
    parameter_range_to_profile_custom,...
    parameter_range_to_profile,...
    data_inclusions]  ...
    = function_load_simulation_settings(simulation_num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WM983b
if	simulation_num == 983101
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 983101 - Scenario 1 - Normoxia
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % For exp data and MLE
    simulation_id = 'Sim983101';
    
    data_file_to_use_array{1} = 'WM983b_confocal_scenario_1_outer.mat';
    data_sheet_to_use_array{1} = '';
    
    p_lower_bounds = [1e-2;1e-2;1e-2;25;25];
    p_upper_bounds = [0.99999;6;1;350;250];
    p_first_guess = [0.9;3;0.5;175;125];
    times_to_use =   [2;3;4;6;8]-2;
    initial_condition = 5000;
    filepath_save = [pwd  '/' simulation_id '/'];
    
    
    % For profile likelihoods
    parameter_range_to_profile_custom = [0;0;0;0;0];
    
    data_inclusions = [1,1,1]; %
    
    
    data_file_to_use_array{2} = 'WM983b_confocal_scenario_1_necrotic.mat';
    data_sheet_to_use_array{2} = '';
    
    data_file_to_use_array{3}= 'WM983b_confocal_scenario_1_inhibited.mat';
    data_sheet_to_use_array{3} = '';

    

elseif simulation_num == 983102
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 983101 - Scenario 2 - Hypoxia
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % For exp data and MLE
    simulation_id = 'Sim983102';
    
    data_file_to_use_array{1} = 'WM983b_confocal_scenario_2_outer.mat';
    data_sheet_to_use_array{1} = '';
    
    p_lower_bounds = [1e-2;1e-2;1e-2;25;25];
    p_upper_bounds = [0.99999;6;1;350;250];
    p_first_guess = [0.9;3;0.5;175;125];
    times_to_use =   [2;4;6;8]-2;
    initial_condition = 5000;
    filepath_save = [pwd  '/' simulation_id '/'];
    
    
    % For profile likelihoods
    parameter_range_to_profile_custom = [0;0;0;0;0];
    
    data_inclusions = [1,1,1]; %
    
    
    data_file_to_use_array{2} = 'WM983b_confocal_scenario_2_necrotic.mat';
    data_sheet_to_use_array{2} = '';
    
    data_file_to_use_array{3}= 'WM983b_confocal_scenario_2_inhibited.mat';
    data_sheet_to_use_array{3} = '';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WM793b
if	simulation_num == 793101
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 793101 - Scenario 1 - Normoxia
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % For exp data and MLE
    simulation_id = 'Sim793101';
    
    data_file_to_use_array{1} = 'WM793b_confocal_scenario_1_outer.mat';
    data_sheet_to_use_array{1} = '';
    
    p_lower_bounds = [1e-2;1e-2;1e-2;25;25];
    p_upper_bounds = [0.99999;6;1;350;250];
    p_first_guess = [0.9;3;0.5;175;125];
    times_to_use =   [2;3;4;6;8]-2;
    initial_condition = 5000;
    filepath_save = [pwd  '/' simulation_id '/'];
    
    
    % For profile likelihoods
    parameter_range_to_profile_custom = [0;0;0;0;0];
    
    data_inclusions = [1,1,1]; %
    
    
    data_file_to_use_array{2} = 'WM793b_confocal_scenario_1_necrotic.mat';
    data_sheet_to_use_array{2} = '';
    
    data_file_to_use_array{3}= 'WM793b_confocal_scenario_1_inhibited.mat';
    data_sheet_to_use_array{3} = '';

    

elseif simulation_num == 793102
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 793101 - Scenario 2 - Hypoxia
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % For exp data and MLE
    simulation_id = 'Sim793102';
    
    data_file_to_use_array{1} = 'WM793b_confocal_scenario_2_outer.mat';
    data_sheet_to_use_array{1} = '';
    
    p_lower_bounds = [1e-2;1e-2;1e-2;25;25];
    p_upper_bounds = [0.99999;6;1;350;250];
    p_first_guess = [0.9;3;0.5;175;125];
    times_to_use =   [2;4;6;8]-2;
    initial_condition = 5000;
    filepath_save = [pwd  '/' simulation_id '/'];
    
    
    % For profile likelihoods
    parameter_range_to_profile_custom = [0;0;0;0;0];
    
    data_inclusions = [1,1,1]; %
    
    
    data_file_to_use_array{2} = 'WM793b_confocal_scenario_2_necrotic.mat';
    data_sheet_to_use_array{2} = '';
    
    data_file_to_use_array{3}= 'WM793b_confocal_scenario_2_inhibited.mat';
    data_sheet_to_use_array{3} = '';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% If would like to calculate the profile likelihood at specific points (not used for manuscript)

if  parameter_range_to_profile_custom(1) == 0 % not using custom then set default values
    parameter_range_to_profile_p1 = (p_lower_bounds(1) + p_upper_bounds(1))./2;
    parameter_range_to_profile_p2 = (p_lower_bounds(2) + p_upper_bounds(2))./2;
    parameter_range_to_profile_p3 = (p_lower_bounds(3) + p_upper_bounds(3))./2;
    parameter_range_to_profile_p4 = (p_lower_bounds(4) + p_upper_bounds(4))./2;
    parameter_range_to_profile_p5 = (p_lower_bounds(5) + p_upper_bounds(5))./2;
    parameter_range_to_profile = {parameter_range_to_profile_p1,parameter_range_to_profile_p2,parameter_range_to_profile_p3,parameter_range_to_profile_p4,parameter_range_to_profile_p5};
end


end