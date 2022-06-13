function [out] = function_Greenspan_2odes(t,y,...
    Rd,...
    lambda,...
    s,...
    Rc,...
    Roinit)

% variables
% y(1) is Outer Radius

% output
out = zeros(1,1);

%% Compute Rn using roots

a3 = 2;
a2 = -3*y(1);
a1=0;
a0 = y(1).*(y(1).^2 - Rc.^2);

polynomial_Rn = [a3 a2 a1 a0];
roots_Rn = roots(polynomial_Rn);

% select root between 0 and Ro, else Rn = 0
Rn_tmp =roots_Rn(logical((roots_Rn<y(1)).*(roots_Rn>0)));

if isempty(Rn_tmp)
    Rn =0;
else
    Rn = Rn_tmp;
end

%% Compute Ri using roots

b3 = -1;
b2 = 0;
b1 = y(1).^2 +2.*(Rn.^3)./y(1) - Rd.^2;
b0 = -2.*Rn.^3; 

polynomial_Ri = [b3 b2 b1 b0];
roots_Ri = roots(polynomial_Ri);

% select root between max(0,Rn) and Ro, else Ri = 0
Ri_tmp =roots_Ri(logical((roots_Ri<y(1)).*(roots_Ri>max(0,Rn))));

if isempty(Ri_tmp)
    Ri =0;
else
    Ri = Ri_tmp;
end


%% dydt(1) - Outer Radius
out(1) = (1./y(1).^2).*(  (s./3).*(y(1).^3 - max(Rn.^3,Ri.^3))   - lambda.*Rn.^3);


end