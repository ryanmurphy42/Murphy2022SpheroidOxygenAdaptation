function out = function_deoxygenation_s(t,s_hypoxia,s_normoxia,tau_s,t_s)

out = s_hypoxia + (s_normoxia - s_hypoxia).*exp(-(1./tau_s).*(t - t_s));

end