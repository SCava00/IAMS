%% Dati iniziali:
clear
clc
% In generale il suffisso _i indica il punto iniziale dove il satellite si
% trova su una singola orbita; _f il punto finale
% Il primo suffisso indica l'orbita, il secondo suffisso indica il punto

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

%% Plot orbita iniziale e finale:
%{
figure(1);
subplot(1,2,1)
plotOrbit(a_i, e_i, i_i, OM_i, om_i)
title('orbita iniziale')

subplot(1,2,2)
plotOrbit(a_f, e_f, i_f, OM_f, om_f);
title('orbita finale')
%}
%% Plot orbite insieme
%{
figure(2);
plotOrbit(a_i, e_i, i_i, OM_i, om_i);
plotOrbit(a_f, e_f, i_f, OM_f, om_f);
title('Orbite')
%}
%% Cambio di piano

% Dati della orbita cambiata di piano:
%%%%%%%%% Suffisso _t1 indica l'orbita cambiata di piano %%%%%%%%%
a_t1 = a_i;
e_t1 = e_i;
i_t1 = i_f;
OM_t1 = OM_f;

[Dv_i_t1,om_t1, th_cp] = changeOrbitalPlane (a_i,e_i,i_i,OM_i,om_i,i_f,OM_f,mu);
th_t1_i = th_cp;
th_i_f = th_cp;
[RR_i_f, ~] = paraorb2rv(a_i, e_i, i_i, OM_i, om_i, th_i_f);

% th_cp = theta del punto di cambio di piano sia sull orbita iniziale che
% sulla prima orbita di trasferiento
% Dv_it1 = costo del cambio di piano

% Calcolo del tempo dal punto iniziale al punto di cambio piano
[Dt_cp] = deltaTime(a_i,e_i, [th_i_i,th_cp]);

%%% Dt_1 = tempo di attesa fino al punto di cambio piano %%%
Dt_1 = Dt_cp;
%%% Dv_1 = costo del cambio piano
Dv_1 = Dv_i_t1;

%% Plot t1 (orbita cambiata di piano) e orbita iniziale
%{
figure(3);
plotOrbit(a_i, e_i, i_i, OM_i, om_i,[th_i_i,th_i_f]);
plotOrbit(a_t1, e_t1, i_t1, OM_t1, om_t1,[th_t1_i,pi]);
title('Orbite cambio piano')
%}

%% Cambio anomalia del pericentro

%%%%%%%%% Suffisso _t2 indica l'orbita cambiata di anomalia del pericentro %%%%%%%%%
% (con anomalia pari a quella dell'orbita finale)

[Dv_2,~,th_rot,Dt_rotom] = orbitalTransfer ('rotom',a_t1,e_t1,[om_t1,om_f]);

%%% Dv_2 = costo del cambio anomalia pericentro

% Calcolo tempi impiegati per passare dal punto iniziale sull'orbita
% cambiata di piano ai punti dove si può far avvenire la manovra di cambio
% anomalia del pericentro:

[Dt_2_caso_a] = deltaTime(a_t1,e_t1, [th_t1_i,th_rot(1,1)]);
[Dt_2_caso_b] = deltaTime(a_t1,e_t1, [th_t1_i,th_rot(2,1)]);


% Assegnamento valori ottenuti:



if Dt_2_caso_a<Dt_2_caso_b
    Dt_2 = Dt_2_caso_a;
    th_t1_f = th_rot(1,1);
    th_t2_i = th_rot(1,2);
else
    Dt_2 = Dt_2_caso_b;
    th_t1_f = th_rot(2,1);
    th_t2_i = th_rot(2,2);
end

%%% Dt_2 = tempo di attesa fino al punto di trasferimento per cambio di om %%%

a_t2 = a_t1;
e_t2 = e_t1;
i_t2 = i_t1;
OM_t2 = OM_t1;
om_t2 = om_f;

%% Plot orbita cambiata di piano (t1) e l'orbita cambiata di om (t2)
%{
figure(4);
plotOrbit(a_t1, e_t1, i_t1, OM_t1, om_t1);
plotOrbit(a_t2, e_t2, i_t2, OM_t2, om_t2);
plotOrbit(a_f, e_f, i_f, OM_f, om_f);
legend('t1','','t2','','Orbita finale')
title('Orbite rifasamento orbitale')
%}

%% Trasferimento bitangente

[DvA,DvB,th_bitan,Dt_bitan,P] = orbitalTransfer ('bitan',[a_t2, a_f], [e_t2,e_f], [om_t2,om_f]);


% *******
% Qui si possono fare due scelte: o tenere il trasferimento [P->A'] che
% costa di meno o quello [A->P'] che ci mette meno tempo (non contando i
% tempi di raggiungimenti dei punti di trasferimento / destinazione)

% *******Utilizzo quello che costa meno solo per comodità:

% Calcolo parametri orbitali dell'orbita di traferimento bitangente:
a_t3 = P(1,1);
e_t3 = P(1,2);
i_t3 = i_t2;
OM_t3 = OM_t2;
om_t3 = om_t2;
th_t3_i = 0;
th_t3_f = pi;

%plotOrbit(a_t3, e_t3, i_t3, OM_t3, om_t3, [th_t3_i,th_t3_f], 1, [0.7176 0.2745 1.0000]);
Dv_3 = DvA(1) + DvA(2);

% Assegnamento valori ottenuti:

th_f_i = th_bitan(2);
th_t2_f = th_bitan(1);

% Calcolo del tempo per raggiungere il punto di manovra
%%% Dt_3 = tempo di attesa fino al punto di trasferimento bitang %%%
[Dt_3] = deltaTime(a_t2,e_t2,[th_t2_i,th_t2_f]);

% Calcolo del tempo di percorrenza dell'orbita di trasferimento
%%% Dt_4 = tempo di trasferimento bitangente %%%
Dt_4 = Dt_bitan(1);

% Calcolo del tempo di percorrenza fino al punto finale
%%% Dt_5 = tempo fino al punto finale %%%
[Dt_5] = deltaTime(a_f,e_f,[th_f_i,th_f_f]);

%% Plot di tutto il trasferimento

% Punto iniziale
plot3(x_i_i,y_i_i,z_i_i,'o','MarkerSize',5, 'LineWidth',2.5,'Color',[0 0.4471 0.7412]);
hold on

% Orbita iniziale:
plotOrbit(a_i, e_i, i_i, OM_i, om_i, [th_i_i,th_i_f], 2, [0 0.4471 0.7412]);

% Punto di cambio piano
plot3(RR_i_f(1),RR_i_f(2),RR_i_f(3),'o','MarkerSize',5, 'LineWidth',2.5,'Color',[0.8510 0.3255 0.0980]);

% Orbita cambiata di piano:
plotOrbit(a_t1, e_t1, i_t1, OM_t1, om_t1, [th_t1_i,th_t1_f], 1, [0.8510 0.3255 0.0980]);

% Orbita riorientata con omega finale
plotOrbit(a_t2, e_t2, i_t2, OM_t2, om_t2, [th_t2_i,th_t2_f], 1, [0.4940 0.1840 0.5560]);

% Orbita di trasferimento bitangente
plotOrbit(a_t3, e_t3, i_t3, OM_t3, om_t3, [th_t3_i,th_t3_f], 1, [0.9290 0.6940 0.1250]);

% Orbita finale
plotOrbit(a_f, e_f, i_f, OM_f, om_f, [th_f_i,th_f_f], 2, [0.6353 0.0784 0.1843]);

% Punto finale
plot3(x_f_f,y_f_f,z_f_f,'o','MarkerSize',5, 'LineWidth',2.5);


legend('Punto iniziale','','','Orbita iniziale','','Punto di cambio piano','Orbita cambiata di piano','','Orbita rifasata con omega finale','','Orbita di trasferimento','','','','Orbita finale','','Punto finale')


%% Costi complessivi

Dt = Dt_1 + Dt_2 + Dt_3 + Dt_4 + Dt_5
Dv = Dv_1 + Dv_2 + Dv_3


