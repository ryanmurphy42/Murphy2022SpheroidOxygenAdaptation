function  [Radii,Ro_plot_times,Rn_plot_times,Ri_plot_times,Rnplus_plot_times] = function_reoxygenation_1simulationbiphase(theta,...
    k,...
    p_infinity,...
    Omega,...
    t_s,...
    plot_times,...
    t_end,...
    variable_id_vec)

% simulate re-oxygenation model

alpha_normoxia = theta(1)/1e7;
alpha_hypoxia = theta(2)/1e7;
tau_alpha  = theta(3);
Rd_normoxia  = theta(4);
Rd_hypoxia  = theta(5);
tau_Rd  = theta(6);
s_normoxia  = theta(7);
s_hypoxia  = theta(8);
tau_s  = theta(9);
lambda_normoxia  = theta(10);
lambda_hypoxia  = theta(11);
tau_lambda  = theta(12);
lambdatilde  = theta(13);
tau_lambdatilde  = theta(14);
outer_radii_at_t0 = theta(15); % at t0
nu = theta(16);

if length(theta) < 17
    necrotic_radii_at_t0 = 0;  % at t0
else
    necrotic_radii_at_t0 = theta(17);  % at t0
end

t_span_hypoxia =  [2,t_s];
t_span_normoxia = [t_s,t_end];

Radii = zeros(length(plot_times),1);
Ro_plot_times = zeros(length(plot_times),1);
Rn_plot_times = zeros(length(plot_times),1);
Ri_plot_times = zeros(length(plot_times),1);
Rnplus_plot_times = zeros(length(plot_times),1);

%% Compute the solution for t<=t_s using Greenspan's model

alpha_at_t = alpha_hypoxia;
p_infinity_2 = 2;

Rc_hypoxia=sqrt((6.*k.*p_infinity_2)./(alpha_at_t.*Omega))./(1e-6);

solution_function_Greenspan_odes = ode15s(@(t,y) function_Greenspan_2odes(t,y,...
    Rd_hypoxia,...
    lambda_hypoxia,...
    s_hypoxia,...
    Rc_hypoxia,...
    outer_radii_at_t0),...
    t_span_hypoxia,...
    outer_radii_at_t0);

%% Compute outer radius and necrotic radius at t_s
y_at_ts = deval(solution_function_Greenspan_odes,t_s);
outer_radii_at_ts_from_greenspan = y_at_ts;

a3 = 2;
a2 = -3*y_at_ts(1);
a1=0;
a0 = y_at_ts(1).*(y_at_ts(1).^2 - Rc_hypoxia.^2);

polynomial_Rn = [a3 a2 a1 a0];
roots_Rn = roots(polynomial_Rn);

% select root between 0 and Ro, else Rn = 0
Rn_tmp =roots_Rn(logical((roots_Rn<y_at_ts(1)).*(roots_Rn>0)));

if isempty(Rn_tmp)
    Rn =0;
else
    Rn = Rn_tmp;
end


necrotic_radii_at_ts_from_greenspan = Rn;

radii_at_ts = [outer_radii_at_ts_from_greenspan % outer radius
    (4*pi/3).*necrotic_radii_at_ts_from_greenspan.^3]; % necrotic core volume

%% Compute the solution for t>t_s using the re-oxygenation model

solution_function_Greenspan_reoxygenation_odes = ode15s(@(t,y) function_reoxygenation_2odes(t,y,...
    alpha_normoxia,...
    alpha_hypoxia,...
    lambda_normoxia,...
    lambda_hypoxia,...
    s_normoxia,...
    s_hypoxia,...
    Rd_normoxia,...
    Rd_hypoxia,...
    tau_alpha,...
    tau_lambda,...
    tau_s,...
    tau_Rd,...
    tau_lambdatilde,...
    lambdatilde,...
    k,...
    p_infinity,...
    Omega,...
    nu,...
    t_s),...
    t_span_normoxia,radii_at_ts);


%% Evaluate the radii at each time step

if max(solution_function_Greenspan_reoxygenation_odes.x) >= max(t_span_normoxia)


    for i=1:length(plot_times)

        t = plot_times(i);

        if t<=t_s

            y=deval(solution_function_Greenspan_odes,t);

            if variable_id_vec(i) == 1
                Ro_plot_times(i)=y(1);
                Radii(i) = y(1);
            end


            if variable_id_vec(i) > 1
                %% Compute Rn using roots

                a3 = 2;
                a2 = -3*y(1);
                a1=0;
                a0 = y(1).*(y(1).^2 - Rc_hypoxia.^2);

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
                    b1 = y(1).^2 +2.*(Rn.^3)./y(1) - Rd_hypoxia.^2;
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


        elseif t> t_s


            alpha_at_t = function_reoxygenation_alpha(t,alpha_hypoxia,alpha_normoxia,tau_alpha,t_s);
            lambda_at_t = function_reoxygenation_lambda(t,lambda_hypoxia,lambda_normoxia,tau_lambda,t_s);
            s_at_t = function_reoxygenation_s(t,s_hypoxia,s_normoxia,tau_s,t_s);
            Rd_at_t = function_reoxygenation_Rd(t,Rd_hypoxia,Rd_normoxia,tau_Rd,t_s);
            lambdatilde_at_t = function_reoxygenation_lambdatilde(t,lambdatilde,tau_lambdatilde,t_s);

            y=deval(solution_function_Greenspan_reoxygenation_odes,t);

            if y(2) < 0
                y(2) = 0;
            end

            if variable_id_vec(i) == 1
                Ro_plot_times(i)=y(1);
                Radii(i)=y(1);
            end

            if variable_id_vec(i) == 2

                if y(2) > 0
                    Rn_plot_times(i)=(3.*y(2)/(4*pi))^(1/3);
                    Radii(i) = (3.*y(2)/(4*pi))^(1/3);
                elseif  y(2)>-1e-6
                    Radii(i) = 0;
                    Rn_plot_times(i)=0;
                end
            end

            if variable_id_vec(i) == 3
                %% Compute Rn+ using roots
                R_c=sqrt((6.*k.*p_infinity)./(alpha_at_t.*Omega))/1e-6;

                a3 = 2;
                a2 = -3*y(1);
                a1=0;
                a0 = y(1).*(y(1).^2 - R_c.^2);

                polynomial_Rnplus = [a3 a2 a1 a0];
                roots_Rnplus = roots(polynomial_Rnplus);

                % select root between 0 and Ro, else Rnplus = 0
                Rnplus_tmp =roots_Rnplus(logical((roots_Rnplus<y(1)).*(roots_Rnplus>0)));

                if isempty(Rnplus_tmp)
                    Rnplus =0;
                else
                    Rnplus = Rnplus_tmp;
                end

                Rnplus_plot_times(i) = Rnplus;

                %% Compute Ri using roots

                b3 = -1;
                b2 = 0;
                b1 = y(1).^2 +2.*(3.*y(2)/(4*pi))./y(1) - Rd_at_t.^2;
                b0 = -2.*(3.*y(2)/(4*pi)); 

                polynomial_Ri = [b3 b2 b1 b0];
                roots_Ri = roots(polynomial_Ri);

                % select root between max(0,Rn) and Ro, else Ri = 0
                Ri_tmp =roots_Ri(logical((roots_Ri<y(1)).*(roots_Ri>max(0,(3.*y(2)/(4*pi))^(1/3)))));

                if isempty(Ri_tmp)
                    Ri =max(0,(3.*y(2)/(4*pi))^(1/3));
                else
                    Ri = Ri_tmp;
                end

                Ri_plot_times(i) = Ri;

                Radii(i) = Ri;
            end

        end

    end

end



end