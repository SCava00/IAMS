function [RR, VV] = paraorb2rv(a, e, i, OM, om, theta)
%
%   Trasformazione da parametri orbitali kepleriani a vettore
%   di stato {R,V} in sistema ECI
%
%   Sistema di riferimento perifocale
%       p = [1;0;0]
%       q = [0;1;0]
%       w = [0;0;1]
%
%   RR = [3x1]
%   VV = [3x1]
%

mu=398600;

% 1) Scrittura di {r,v} in sistema perifocale

p=a*(1-e^2);                % Calcolo semilato retto

r=p/(1+e*cos(theta));       % Calcolo modulo del raggio

rr(1)=r*cos(theta);         % Calcolo componenti di r in PF
rr(2)=r*sin(theta);
rr(3)=0;
rr=rr';

vv(1)=-sqrt(mu/p)*sin(theta);       % Calcolo componenti di v in PF
vv(2)=sqrt(mu/p)*(e+cos(theta));
vv(3)=0;
vv=vv';

% 2) Costruzione matrici di rotazione

R_OM=[cos(OM) sin(OM) 0;
      -sin(OM) cos(OM) 0;
      0 0 1];

R_om=[cos(om) sin(om) 0;
      -sin(om) cos(om) 0;
      0 0 1];

R_i=[1 0 0;
    0 cos(i) sin(i);
    0 -sin(i) cos(i)];

R_pf2eci=R_OM'*(R_i'*R_om');

% 3) Calcolo di R e V in ECI

RR=R_pf2eci*rr;
VV=R_pf2eci*vv;

