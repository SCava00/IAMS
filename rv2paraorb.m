function [a, e, i, OM, om, th] = rv2paraorb(RR, VV, mu)
%
%   Trasformazione da vettore di stato {R,V} in sistema ECI, a parametri
%   kepleriani
%
%   RR = [3x1]
%   VV = [3x1]
%   mu = [1x1]
%
%   Sistema di riferimento ECI
%       I = [1;0;0]
%       J = [0;1;0]
%       K = [0;0;1]
%

% 1) Calcolo norme vettori

R=norm(RR);
V=norm(VV);

% 2) Calcolo di a

E=(V^2)/2-mu/R;
a=-mu/(2*E);

% 3) Calcolo e

hh=cross(RR,VV);
h=norm(hh);

ee=(cross(VV,hh))./mu - RR./R;
e=norm(ee);

% 4) Calcolo i

i=acos(hh(3)/h);

% 5) Calcolo OM

K=[0;0;1];
NN=(cross(K,hh))./(norm(cross(K,hh)));

if(NN(2)>0)
    OM=acos(NN(1));
end

if(NN(2)<0)
    OM=2*pi-acos(NN(1));
end

% 6) Calcolo om

if(ee(3)>0)
    om=acos(dot(ee,NN)/e);
end
if(ee(3)<0)
    om=2*pi-acos(dot(ee,NN)/e);
end

% 7) Calcolo theta

Vr=dot(RR,VV)/R;
if(Vr>0)
    th=acos(dot(RR,ee)/(R*e));
end
if(Vr<0)
    th=2*pi-acos(dot(RR,ee)/(R*e));
end


