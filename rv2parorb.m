function [a, e, i, OM, om, th] = rv2parorb(rr, vv, mu)
%  Conversione da vettori di stato ([x, y, z], [vx, vy, vz]) in forma vettoriale a parametri kepleriani (a, e i, ohm, omega, theta)
% Parametro extra di input che definisce la costante di gravitazionale

r = norm(rr);   %Modulo del raggio
v = norm(vv);   %Modulo del vettore velocità

a = mu/( (2*mu/r) - v^2 );  %Calcolo semiasse maggiore tramite la formula dell'energia dell'orbita

hh = cross(rr, vv); %Calcolo il vettore momento angolare specifico

h = norm(hh);       %Modulo vettore momento angolare

ee = (cross(vv, hh)/mu - (rr/r)); %Calcolo vettore eccentricità
e = norm(ee);


i = acos(hh(3)/h); %Inclinazione piano orbitale [rad]

nn = cross([0 0 1]', hh)./norm(cross([0 0 1]', hh)); %Linea dei nodi


if nn(2)>= 0
    OM = acos(nn(1));
else
    OM = 2*pi() - acos(nn(1));
end
% Calcolo dell'ascensione retta del nodo ascendente


if ee(3) >= 0
    om = acos(ee'*nn/e);
else
    om = 2*pi() - acos(ee'*nn/e);
end

% Calcolo argomento del pericentro

vr = dot(vv, rr)/r;

if vr >= 0
    th = acos(dot(ee, rr)/(norm(ee)*norm(rr)));
else
    th = 2*pi() - acos(dot(ee, rr)/(norm(ee)*norm(rr)));
end

end
