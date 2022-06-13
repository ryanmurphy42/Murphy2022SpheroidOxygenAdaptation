function out = function_deoxygenation_lambda(t,lambda_hypoxia,lambda_normoxia,tau_lambda,t_s)

out = lambda_hypoxia + (lambda_normoxia - lambda_hypoxia).*exp(-(1./tau_lambda).*(t - t_s));

end