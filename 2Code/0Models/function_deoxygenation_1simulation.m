function  [Radii,Ro_plot_times,Rn_plot_times,Ri_plot_times,Rnplus_plot_times] = function_deoxygenation_1simulation(theta,...
    k,...
    p_infinity,...
    Omega,...
    t_s,...
    plot_times,...
    t_end,...
    variable_id_vec)

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
lambdahat  = theta(13);
tau_lambdahat  = theta(14);
outer_radii_at_ts = theta(15);
necrotic_radii_at_ts = 0;

radii_at_ts = [outer_radii_at_ts
    necrotic_radii_at_ts];

t_span_hypoxia =  [t_s,t_end];

Radii = zeros(length(plot_times),1);
Ro_plot_times = zeros(length(plot_times),1);
Rn_plot_times = zeros(length(plot_times),1);
Ri_plot_times = zeros(length(plot_times),1);
Rnplus_plot_times = zeros(length(plot_times),1);


solution_function_deoxygenation_2odes = ode15s(@(t,y) function_deoxygenation_2odes(t,y,...
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
    tau_lambdahat,...
    lambdahat,...
    k,...
    p_infinity,...
    Omega,...
    t_s),...
    t_span_hypoxia,radii_at_ts);

if max(solution_function_deoxygenation_2odes.x) >= max(t_span_hypoxia)


    for i=1:length(plot_times)

        t = plot_times(i);
        alpha_at_t = function_deoxygenation_alpha(t,alpha_hypoxia,alpha_normoxia,tau_alpha,t_s);
        lambda_at_t = function_deoxygenation_lambda(t,lambda_hypoxia,lambda_normoxia,tau_lambda,t_s);
        s_at_t = function_deoxygenation_s(t,s_hypoxia,s_normoxia,tau_s,t_s);
        Rd_at_t = function_deoxygenation_Rd(t,Rd_hypoxia,Rd_normoxia,tau_Rd,t_s);
        lambdahat_at_t = function_deoxygenation_lambdahat(t,lambdahat,tau_lambdahat,t_s);

        y=deval(solution_function_deoxygenation_2odes,t);

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

else

end


end