function [rr, vv] = parorb2rv(a, e, i, ohm, omega, theta, mu)
%
%
Ro = [cos(ohm) sin(ohm) 0
     -sin(ohm) cos(ohm) 0
     0 0 1];
 
Ri = [1 0 0
     0 cos(i) sin(i)
     0 -sin(i) cos(i)];
 
Romega = [cos(omega) sin(omega) 0
         -sin(omega) cos(omega) 0
         0 0 1];

Tpfge = (Romega*Ri*Ro)';


p = a*(1-e^2);

r = p / (1 + e*cos(theta));
rpf = [r*cos(theta) r*sin(theta) 0]';

vpf = [-sqrt(mu/p)*sin(theta) sqrt(mu/p)*(e+cos(theta)) 0]';

rr = Tpfge * rpf;
vv = Tpfge * vpf;

end

