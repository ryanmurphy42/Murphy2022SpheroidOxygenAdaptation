function out = function_deoxygenation_alpha(t,alpha_hypoxia,alpha_normoxia,tau_alpha,t_s)

out = alpha_hypoxia + (alpha_normoxia - alpha_hypoxia).*exp(-(1./tau_alpha).*(t - t_s));

end