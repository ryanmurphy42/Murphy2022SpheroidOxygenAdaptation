function out = function_reoxygenation_Rd(t,Rd_hypoxia,Rd_normoxia,tau_Rd,t_s)

out = Rd_normoxia  + (Rd_hypoxia - Rd_normoxia).*exp(-(1./tau_Rd).*(t - t_s));

end