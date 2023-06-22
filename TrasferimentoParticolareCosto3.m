%% Dati iniziali:
clear
clc
% In generale il suffisso _i indica il punto iniziale dove il satellite si
% trova su una singola orbita; _f il punto finale
% Il primo suffisso indica l'orbita, il secondo suffisso indica il punto
% Il primo suffisso indica l'orbita e il secondo indica se si riferisce a
% punto iniziale o finale

% Punto iniziale:
%%%%%%%%% Suffisso _i indica l'orbita iniziale %%%%%%%%%
x_i_i = -5919.8013;
y_i_i = -24.7819; 
z_i_i = 5419.0110;
vx_i_i = -2.9010;
vy_i_i = -6.0330;
vz_i_i = -2.3130;

% Punto finale:
%%%%%%%%% Suffisso _f indica l'orbita finale %%%%%%%%%
a_f = 13950.0000; 
e_f = 0.2715;
i_f = 0.9413;
OM_f = 1.1850;
om_f = 1.3280;
th_f_f = 1.8460;

mu = 398600;

%% Calcolo parametri orbitali dell'orbita iniziale:

RR = [x_i_i; y_i_i; z_i_i];
VV = [vx_i_i; vy_i_i; vz_i_i];
[a_i, e_i, i_i, OM_i, om_i, th_i_i] = rv2paraorb(RR, VV, mu);

%% Calcolo parametri spaziali punto finale

[RR, ~]  = paraorb2rv(a_f, e_f, i_f, OM_f, om_f, th_f_f);

% Assegnamento valori ottenuti 
x_f_f = RR(1);
y_f_f = RR(2);
z_f_f = RR(3);

%% Descrizione trasferimento:

% Ottimizzazione costo
% 1) Orbita iniziale
% 2) Aspettare fino al punto 'simmetrico' rispetto al punto di intersezione
%    tra l'orbita finale e il piano iniziale più lontano
% 3) Iniziare un trasferimento 'alla hohmann' fino al punto trovato, dove
%    l'orbita finale è un orbita circolare con raggio uguale al raggio del
%    punto individuato e l'orbita iniziale è un orbita circolare di raggio
%    pari alla distanza del punto dove far avvenire la manovra (per fare i
%    calcoli utilizzare i delta delle velocità assolute)
% 4) Far avvenire il cambio di piano non standard
% 5) Aspettare fino al punto finale

%% 1) Orbita iniziale
Dv_1 = 0;
Dt_1 = 0;

% Punto iniziale
plot3(x_i_i,y_i_i,z_i_i,'o','MarkerSize',5, 'LineWidth',2.5);
hold on;

%% 2) Primo punto di manovra
% per trovare il punto di intersezione delle due orbite utilizzo
% changeOrbital plane per trovare th dell'orbita finale
[~,om_teorico,th_f_i] = changeOrbitalPlane (a_f,e_f,i_f,OM_f,om_f,i_i,OM_i,mu);

% ora per ricavare la theta al quale far avvenire la manovra devo
% considerare il DOM che avrebbe il cambio di piano
Dth = om_teorico - om_f;
th_i_f = (th_f_i + Dth) - pi;
% ma questo sarebbe il theta se l'orbita iniziale e finale avessero la
% stessa om; quindi rispetto all'orbita iniziale il th sarà:
Dom_i_f = om_f - om_i;
th_i_f = th_i_f + Dom_i_f;

plotOrbit(a_i,e_i,i_i,OM_i,om_i,[th_i_i,th_i_f],2,[0 0.4471 0.7412]);


Dv_2 = 0;
[Dt_2] = deltaTime(a_i,e_i,[th_i_i,th_i_f]);

%% 3) Trasferimento alla hohmann
% ora ricavo le coordinate spaziali del punto iniziale
[RR, VV] = paraorb2rv(a_i, e_i, i_i, OM_i, om_i, th_i_f);
plot3(RR(1),RR(2),RR(3),'o','MarkerSize',5, 'LineWidth',2.5)
% ricavo così il raggio della circonferenza da cui fare il trasferimento
r_i = sqrt(RR(1)^2 + RR(2)^2 + RR(3)^2);
% e ricavo anche la velocità in coordinate
vx_i_f = VV(1);
vy_i_f = VV(2);
vz_i_f = VV(3);

% ora ricavo le coordinate spaziali del punto finale
[RR, ~] = paraorb2rv(a_f, e_f, i_f, OM_f, om_f, th_f_i);
plot3(RR(1),RR(2),RR(3),'o','MarkerSize',5, 'LineWidth',2.5)
% ricavo così il raggio della circonferenza a cui fare il trasferimento
r_f = sqrt(RR(1)^2 + RR(2)^2 + RR(3)^2);

[DvA,DvB,th,Dt,P] = orbitalTransfer ('hohmann',[r_i,r_f],0,0,0);
a_t1 = P(1);
e_t1 = P(2);
i_t1 = i_i;
OM_t1 = OM_i;
om_t1 = om_i +th_i_f;
th_t1_i = 0;
th_t1_f = pi;

plotOrbit(a_t1,e_t1,i_t1,OM_t1,om_t1,[th_t1_i,th_t1_f],1,[0.9290 0.6940 0.1250]);

% per il calcolo del costo faccio la differenza tra la velocità che si ha
% sull'orbita iniziale e quella sulla di trasferimento alla hohmann
[~, VV] = paraorb2rv(a_t1, e_t1, i_t1, OM_t1, om_t1, th_t1_i);
vx_t1_i = VV(1);
vy_t1_i = VV(2);
vz_t1_i = VV(3);
Dvx = vx_i_f - vx_t1_i;
Dvy = vy_i_f - vy_t1_i;
Dvz = vz_i_f - vz_t1_i;

Dv_3 = sqrt(Dvx^2 + Dvy^2 + Dvz^2);
Dt_3 = Dt;

%% 6) Cambio di piano
% Calcolo le velocità nel punto di manovra sulle due orbite
[~, VV] = paraorb2rv(a_t1, e_t1, i_t1, OM_t1, om_t1, th_t1_f);
vx_t1_f = VV (1);
vy_t1_f = VV (2);
vz_t1_f = VV (3);

[~, VV] = paraorb2rv(a_f, e_f, i_f, OM_f, om_f, th_f_i);
vx_f_i = VV (1);
vy_f_i = VV (2);
vz_f_i = VV (3);

Dvx = vx_f_i - vx_t1_f;
Dvy = vy_f_i - vy_t1_f;
Dvz = vz_f_i - vz_t1_f;

Dv_4 = sqrt(Dvx^2 + Dvy^2 + Dvz^2);
Dt_4 = 0;

%% 5) Raggiungimento punto finale

Dv_5 = 0;
[Dt_5] = deltaTime(a_f,e_f,[th_f_i,th_f_f]);

% Punto finale
plot3(x_f_f,y_f_f,z_f_f,'o','MarkerSize',5, 'LineWidth',2.5);

%% Costi
Dt = Dt_1 + Dt_2 + Dt_3 + Dt_4 + Dt_5
Dv = Dv_1 + Dv_2 + Dv_3 + Dv_4 + Dv_5

plotOrbit(a_f,e_f,i_f,OM_f,om_f,[th_f_i,th_f_f],2,[0.6353 0.0784 0.1843]);

legend('Punto iniziale','','','Orbita iniziale','','Punto prima manovra','Punto di cambio piano','Orbita di trasferimento','','Punto finale','','','Orbita finale')

% NOTA: si nota dal plot che la il punto iniziale è gia vicinissimo alla
% orbita di trasferimento quindi volendo si può pensare di progettare un
% trasferimento che faccia risparmiare tempo 
