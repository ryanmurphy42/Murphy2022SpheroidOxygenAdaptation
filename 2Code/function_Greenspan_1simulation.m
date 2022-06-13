function [Radii,Ro_plot_times,Rn_plot_times,Ri_plot_times] = function_Greenspan_1simulation(parameters,...
    plot_times,...
    variable_id_vec,...
    t_end)

% variable_id_vec = 1 - outer, 2- necrotic, 3 - inhibited

Rd = parameters(1);
lambda = parameters(2);
s = parameters(3);
Rc = parameters(4);
Roinit = parameters(5);

% simulate Greenspans model

t_span =  [0,t_end];

Radii = zeros(length(plot_times),1);


solution_function_Greenspan_2odes = ode15s(@(t,y) function_Greenspan_2odes(t,y,...
    Rd,...
    lambda,...
    s,...
    Rc,...
    Roinit),t_span,Roinit);

for i=1:length(plot_times)


    t = plot_times(i);
    y=deval(solution_function_Greenspan_2odes,t);

    if variable_id_vec(i) == 1
        Ro_plot_times(i)=y(1);
        Radii(i) = y(1);

    end

    if variable_id_vec(i) > 1
        %% Compute Rn using roots

        a3 = 2;
        a2 = -3*y(1);
        a1=0;
        a0 = y(1).*(y(1).^2 - Rc.^2);

        polynomial_Rn = [a3 a2 a1 a0];
        roots_Rn = roots(polynomial_Rn);

        % select root between 0 and Ro, else Rn = 0
        Rn_tmp =roots_Rn(logical((roots_Rn<y(1)).*(roots_Rn>0)));

        if isempty(Rn_tmp)
            Rn =0;
        else
            Rn = Rn_tmp;
        end

        if variable_id_vec(i) == 2
            Rn_plot_times(i) = Rn;
            Radii(i)=Rn;
        end

        if variable_id_vec(i)  == 3
            %% Compute Ri using roots

            b3 = -1;
            b2 = 0;
            b1 = y(1).^2 +2.*(Rn.^3)./y(1) - Rd.^2;
            b0 = -2.*Rn.^3;

            polynomial_Ri = [b3 b2 b1 b0];
            roots_Ri = roots(polynomial_Ri);

            % select root between max(0,Rn) and Ro, else Ri = 0
            Ri_tmp =roots_Ri(logical((roots_Ri<y(1)).*(roots_Ri>max(0,Rn))));

            if isempty(Ri_tmp)
                Ri = 0;
            else
                Ri = Ri_tmp;
            end

            Ri_plot_times(i) = Ri;
            Radii(i) = Ri;
        end
    end

end


end