function [out] = function_deoxygenation_2odes(t,y,...
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
    t_s)

% variables
% y(1) is Outer Radius
% y(2) is Necrotic VOLUME

% output
out = zeros(2,1);

%% Compute alpha, lambda, s, and w at time t
alpha_at_t = function_deoxygenation_alpha(t,alpha_hypoxia,alpha_normoxia,tau_alpha,t_s);
lambda_at_t = function_deoxygenation_lambda(t,lambda_hypoxia,lambda_normoxia,tau_lambda,t_s);
s_at_t = function_deoxygenation_s(t,s_hypoxia,s_normoxia,tau_s,t_s);
Rd_at_t = function_deoxygenation_Rd(t,Rd_hypoxia,Rd_normoxia,tau_Rd,t_s);
lambdahat_at_t = function_deoxygenation_lambdahat(t,lambdahat,tau_lambdahat,t_s);

%% Compute Rn+ using roots

% y(1)
% y(2)

R_c=sqrt((6.*k.*p_infinity)./(alpha_at_t.*Omega))./(1e-6);

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
    Ri =0;
else
    Ri = Ri_tmp;
end


%% dydt(1) - Outer Radius
out(1) = (1./y(1).^2).*(  (s_at_t./3).*(y(1).^3 - max((3.*y(2)/(4*pi)),Ri.^3))   - lambda_at_t.*(3.*y(2)/(4*pi)));

%% dydt(2) - Necrotic VOLUME
out(2) = 3.*lambdahat_at_t.*( (4.*pi./3).*Rnplus.^3 - y(2) ) -3.*lambda_at_t.*y(2);



end