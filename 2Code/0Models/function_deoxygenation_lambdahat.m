function out = function_deoxygenation_lambdahat(t,lambdahat,tau_lambdahat,t_s)

lambdahat_max = 20000000;

if (lambdahat.*exp((1./tau_lambdahat).*(t-t_s))) > lambdahat_max
    out = lambdahat_max;
else
    out = lambdahat.*exp((1./tau_lambdahat).*(t-t_s));
end

end