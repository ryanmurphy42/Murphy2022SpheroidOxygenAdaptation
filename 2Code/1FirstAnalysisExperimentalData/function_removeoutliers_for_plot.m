function [scenario_1_days_to_plot_outer_tmp2,scenario_1_plotdata_outer_tmp2] = function_removeoutliers_for_plot(scenario_1_days_to_plot_outer_tmp,scenario_1_plotdata_outer_tmp)

unique_scenario_1_days_to_plot_outer_tmp = sort(unique(scenario_1_days_to_plot_outer_tmp));
for i=1:length(unique_scenario_1_days_to_plot_outer_tmp)
    outliers_scenario_1_plotdata_outer_tmp = isoutlier(scenario_1_plotdata_outer_tmp(scenario_1_days_to_plot_outer_tmp==unique_scenario_1_days_to_plot_outer_tmp(i)),'quartiles');
    sum_outliers_scenario_1_plotdata_outer_tmp = sum(outliers_scenario_1_plotdata_outer_tmp);
    if sum_outliers_scenario_1_plotdata_outer_tmp>0
        %if there are outliers then remove them from the plot
        scenario_1_plotdata_outer_tmp2 = scenario_1_plotdata_outer_tmp(scenario_1_days_to_plot_outer_tmp==unique_scenario_1_days_to_plot_outer_tmp(i));
        for j=1:length(outliers_scenario_1_plotdata_outer_tmp)
            if outliers_scenario_1_plotdata_outer_tmp(j) > 0 %if an outlier
                % find radius value
                scenario_1_plotdata_outer_value =  scenario_1_plotdata_outer_tmp2(j);
                % find in data
                while sum(      (scenario_1_days_to_plot_outer_tmp  == unique_scenario_1_days_to_plot_outer_tmp(i)).*(scenario_1_plotdata_outer_tmp  == scenario_1_plotdata_outer_value) ) >0
                    % remove from tmp
                    [find_1,~] = find((scenario_1_days_to_plot_outer_tmp  == unique_scenario_1_days_to_plot_outer_tmp(i)).*(scenario_1_plotdata_outer_tmp  == scenario_1_plotdata_outer_value));
                    if find_1 == 1
                        scenario_1_days_to_plot_outer_tmp = [scenario_1_days_to_plot_outer_tmp(2:end)];
                        scenario_1_plotdata_outer_tmp = [scenario_1_plotdata_outer_tmp(2:end)];
                    elseif find_1 == length(scenario_1_days_to_plot_outer_tmp)
                        scenario_1_days_to_plot_outer_tmp = [scenario_1_days_to_plot_outer_tmp(1:end-1)];
                        scenario_1_plotdata_outer_tmp = [scenario_1_plotdata_outer_tmp(1:end-1)];
                    else
                        scenario_1_days_to_plot_outer_tmp = [scenario_1_days_to_plot_outer_tmp(1:find_1-1);scenario_1_days_to_plot_outer_tmp(find_1+1:end)];
                        scenario_1_plotdata_outer_tmp = [scenario_1_plotdata_outer_tmp(1:find_1-1);scenario_1_plotdata_outer_tmp(find_1+1:end)];
                    end
                end
            end
        end
    end
end

scenario_1_days_to_plot_outer_tmp2=scenario_1_days_to_plot_outer_tmp;
scenario_1_plotdata_outer_tmp2 = scenario_1_plotdata_outer_tmp;

end