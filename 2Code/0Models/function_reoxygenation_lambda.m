function out = function_reoxygenation_lambda(t,lambda_hypoxia,lambda_normoxia,tau_lambda,t_s)

out = lambda_normoxia + (lambda_hypoxia - lambda_normoxia).*exp(-(1./tau_lambda).*(t - t_s));

end