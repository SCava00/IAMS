function [Dv] = deltaVtang (a_i,a_f,r)
%
%       Funzione che dato un raggio r e il semiasse maggiore di un orbita
%       iniziale e una finale, restituisce la differenza di velocità da
%       applicare.
%       Funziona solo nei punti apsidali, perchè la velocità è nella stessa
%       direzione
%

mu=398600;

v_i=sqrt(2*mu* ( (1/r) - 1/(2*a_i) ) );
v_f=sqrt(2*mu* ( (1/r) - 1/(2*a_f) ) );
Dv=abs(v_f-v_i);