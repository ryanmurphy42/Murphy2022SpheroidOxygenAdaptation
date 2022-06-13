function out = function_reoxygenation_alpha(t,alpha_hypoxia,alpha_normoxia,tau_alpha,t_s)

out = alpha_normoxia  + (alpha_hypoxia - alpha_normoxia).*exp(-(1./tau_alpha).*(t - t_s));

end