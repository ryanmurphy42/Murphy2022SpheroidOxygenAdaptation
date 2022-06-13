function out = function_reoxygenation_s(t,s_hypoxia,s_normoxia,tau_s,t_s)

out = s_normoxia  + (s_hypoxia - s_normoxia).*exp(-(1./tau_s).*(t - t_s));

end