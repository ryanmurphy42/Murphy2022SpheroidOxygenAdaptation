function out = function_deoxygenation_Rd(t,Rd_hypoxia,Rd_normoxia,tau_Rd,t_s)

out = Rd_hypoxia + (Rd_normoxia - Rd_hypoxia).*exp(-(1./tau_Rd).*(t - t_s));

end