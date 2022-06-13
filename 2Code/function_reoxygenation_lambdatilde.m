function out = function_reoxygenation_lambdatilde(t,lambdatilde,tau_lambdatilde,t_s)

lambdatilde_max = 20000000;

if (lambdatilde.*exp((1./tau_lambdatilde).*(t-t_s))) > lambdatilde_max
    out = lambdatilde_max;
else
    out = lambdatilde.*exp((1./tau_lambdatilde).*(t-t_s));
end

end