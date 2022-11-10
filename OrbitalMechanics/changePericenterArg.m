function [Dv, thi, thf] = changePericenterArg(a, e, omi, omf)
%CHANGEPERICENTERARG
%   

mu = 398600;

Dom = abs(omf-omi);

thi = Dom/2;
thf =  - thi;

Dv = 2 * sqrt(mu/(a*(1-e^2)))*e*sin(Dom/2);

end

