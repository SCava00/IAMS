function [rr, vv] = parorb2rv(a, e, i, OM, om, th, mu)
%
%
Ro = [cos(OM) sin(OM) 0
     -sin(OM) cos(OM) 0
     0 0 1];
 
Ri = [1 0 0
     0 cos(i) sin(i)
     0 -sin(i) cos(i)];
 
Romega = [cos(om) sin(om) 0
         -sin(om) cos(om) 0
         0 0 1];

Tpfge = (Romega*Ri*Ro)';


p = a*(1-e^2);

r = p / (1 + e*cos(th));
rpf = [r*cos(th) r*sin(th) 0]';

vpf = [-sqrt(mu/p)*sin(th) sqrt(mu/p)*(e+cos(th)) 0]';

rr = Tpfge * rpf;
vv = Tpfge * vpf;

end

