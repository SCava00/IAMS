function plotOrbit(a, e, i, OM, om, Dth)
%
%   Funzione che disegna l'orbita assegnata nel sistema ECI
%   L'origine del grafico è fissata nel fuoco dell'orbita
%
%   Dth è un vettore 1x2 contenente la theta iniziale e finale
%   dell'orbita da disegnare
%
%   La sfera è la terra
%
%

%{
a=20000;
e=0.5;
i=pi/4;
om=pi/4;
OM=0;
Dth=[pi,pi/2];
%}

if(nargin==5)
    Dth=[0,2*pi];   % Se Dth non immesso fa il plot di tutta l'ellise
end

prec=50;

% Controllo che theta finale sia maggiore di theta iniziale e eventuale
% spostamento di Dth(1) di 2pi
while (Dth(2)<Dth(1))
    Dth(1)=Dth(1)-2*pi;
end

R_pian=6371;                             % Raggio pianeta (terra)
theta=linspace(Dth(1),Dth(2),prec);      % Vettore theta che va da tho a thf

% Creazione vettori dell'orbita
Rx=[zeros(prec,1)];
Ry=[zeros(prec,1)];
Rz=[zeros(prec,1)];
val=1;
for th=theta(1:prec)      % Calcolo di R in sistema ECI
    [R, ~] = paraorb2rv(a, e, i, OM, om, th); 
    Rx(val)=[R(1)];       % Creazione di vettori di coordinate
    Ry(val)=[R(2)];
    Rz(val)=[R(3)];
    val=val+1;
end



% Plot della proiezione dell orbita sulla terra
PRx=[zeros(prec,1)];
PRy=[zeros(prec,1)];
PRz=[zeros(prec,1)];
var=1;
for th=theta(1:prec)      
    [R, ~] = paraorb2rv(R_pian, 0, i, OM, om, th); 
    PRx(var)=[R(1)];       
    PRy(var)=[R(2)];
    PRz(var)=[R(3)];
    var=var+1;
end


%{
% Parametrizzazione sfera
[Xs,Ys,Zs]=sphere;
Xs=Xs*R_pian;
Ys=Ys*R_pian;
Zs=Zs*R_pian;
%}


plot3(Rx,Ry,Rz,'LineWidth',1.5)     % Plot orbita
hold on;
axis equal;
plot3(PRx,PRy,PRz,'LineWidth',1.5)  % Plot proiezione orbita su terra
earth_sphere;           % Funzione che plotta la terra
%surf(Xs,Ys,Zs,'FaceAlpha',0.5)      % Plot sfera

