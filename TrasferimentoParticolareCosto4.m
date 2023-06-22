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
% 1) Orbita iniziale, punto iniziale
% 2) Subito do un Delta_V per ampliare l'orbita fino ad arrivare a contatto
%    con l'orbita che trasferisce al punto dove farò avvenire il cambio di
%    piano (contemporanteo con l'ingresso sull'orbita finale)
% 3) Trasferirsi con cambio di orbita istantaneo sulla orbita di
%    trasferimento 'alla hohmann' che porta al punto di cambio piano più
%    lontano
% 4) Far avvenire il cambio di piano non standard
% 5) Aspettare fino al punto finale

%% 1) Orbita iniziale

% Punto iniziale
plot3(x_i_i,y_i_i,z_i_i,'o','MarkerSize',5, 'LineWidth',2.5);
hold on;
plotOrbit(a_i,e_i,i_i,OM_i,om_i,[0,2*pi],3, [0 0.4471 0.7412]);


%% 1 / 2) Impulso iniziale e ingresso su orbita hohmann
% Scrivo i valori dell'orbita di trasferimento alla hohmann già trovati con
% il caso precedente
a_t2 = 1.269844438942552e+04;
e_t2 = 0.396275203388309;
i_t2 = 0.890171516561191;
OM_t2 = 0.839030517654833;
om_t2 = 1.494455299846049;
th_t2_f = pi; 


% Creo un versore del vettore posizione iniziale del satellite, in questo
% modo posso applicare un Dv come Dv = n * versore; dove n posso sceglierlo
% ed ottimizzarlo in modo che sia il minore possibile che permetta di
% raggiungere l'orbita di trasferimento alla hohmann

RR = [x_i_i; y_i_i; z_i_i];
VV = [vx_i_i; vy_i_i; vz_i_i];
vers_Ri = RR./norm(RR);
vers_Vi = VV./norm(VV);

% Scelta di n e quindi del quantitativo di Dv iniziale
n = 0.8058;
Dv_1 = n*(vers_Ri+vers_Vi);

% Ridefinisco quindi l'orbita che avrà inizialmente il satellite calcolando
% i nuovi parametri orbitali
VV_nuova = VV + Dv_1;
Dv_1 = norm(Dv_1);

[a_t1,e_t1,i_t1,OM_t1,om_t1,th_t1_i] = rv2paraorb(RR, VV_nuova, mu);


% Devo quindi trovare i punti di intersezione tra la nuova orbita iniziale
% e quella del trasferimento alla hohmann

Dom = -om_t1+om_t2; % definisco l'angolo di differenza tra le due orbite
prec = 100000; % imposto la precisione
theta = linspace(0,pi,prec);
R1 = zeros(1,prec);
R2 = zeros(1,prec);
k = 1;
inters = zeros(1,2); % prealloco i valori di intersezione
for i = 1:prec
    [RR, ~] = paraorb2rv(a_t1,e_t1,i_t1,OM_t1,om_t1,theta(i)+Dom);
    R1(i) = norm(RR);
    [RR, ~] = paraorb2rv(a_t2,e_t2,i_t2,OM_t2,om_t2,theta(i));
    R2(i) = norm(RR);
    if (i>1 && ((R1(i-1)>R2(i-1) && R1(i)<=R2(i)) || (R1(i-1)<R2(i-1) && R1(i)>=R2(i) ) ))
        inters(k) = theta(i); % punti di intersezione delle due orbite (theta rispetto a t2)
        k = k+1;
    end
end

% ora devo calcolare il costo per ognuno dei due punti di intersezione

% costo nel primo punto:
[~, VV1] = paraorb2rv(a_t1,e_t1,i_t1,OM_t1,om_t1,inters(1)+Dom);
[~, VV2] = paraorb2rv(a_t2,e_t2,i_t2,OM_t2,om_t2,inters(1));
Dv_a = norm(VV1-VV2);

% costo nel secondo punto:
[~, VV1] = paraorb2rv(a_t1,e_t1,i_t1,OM_t1,om_t1,inters(2)+Dom);
[~, VV2] = paraorb2rv(a_t2,e_t2,i_t2,OM_t2,om_t2,inters(2));
Dv_b = norm(VV1-VV2);

% scelgo quindi il caso che mi permette il costo minore
if (Dv_a<Dv_b)
    th_t1_f = Dom + inters(1);
    th_t2_i = inters(1);
    Dv_2 = Dv_a;
else
    th_t1_f = Dom + inters(2);
    th_t2_i = inters(2);
    Dv_2 = Dv_b;
end



% plot delle traiettorie:
plotOrbit(a_t1, e_t1, i_t1, OM_t1, om_t1,[th_t1_i th_t1_f],1,[0.8510 0.3255 0.0980]);
plotOrbit(a_t2,e_t2,i_t2,OM_t2,om_t2,[th_t2_i th_t2_f],1,[0.9290 0.6940 0.1250]);

% Calcolo del tempo di permanenza su t1
[Dt_1] = deltaTime(a_t1,e_t1,[th_t1_i,th_t1_f]);

% Calcolo del tempo di permanenza su t2
[Dt_2] = deltaTime(a_t2,e_t2,[th_t2_i,th_t2_f]);

%{
figure
hold off
plot(theta,R1);
hold on;
plot(theta,R2);
%}

%% 3) Cambio di piano
th_f_i = 3.095882922908429; % calcolato in 'TrasferimentoparticolareCosto3'

% Calcolo le velocità nel punto di manovra sulle due orbite
[~, VV1] = paraorb2rv(a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2_f);

[~, VV2] = paraorb2rv(a_f, e_f, i_f, OM_f, om_f, th_f_i);

Dv_3 = norm(VV1-VV2);
Dt_3 = 0;

%% 4) Raggiungimento punto finale

Dv_4 = 0;
[Dt_4] = deltaTime(a_f,e_f,[th_f_i,th_f_f]);


%% Costi
Dt = Dt_1 + Dt_2 + Dt_3 + Dt_4
Dv = Dv_1 + Dv_2 + Dv_3 + Dv_4


plotOrbit(a_f,e_f,i_f,OM_f,om_f,[th_f_i,th_f_f],2,[0.6353 0.0784 0.1843]);

% Punto finale
plot3(x_f_f,y_f_f,z_f_f,'o','MarkerSize',5, 'LineWidth',2.5,'Color',[0.6353 0.0784 0.1843]);

legend('Punto iniziale','Orbita iniziale','','Prima orbita di trasferimento','','Seconda orbita di trasferimento','','','','Orbita finale','','Punto finale')

