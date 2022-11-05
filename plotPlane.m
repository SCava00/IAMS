function plotPlane(a,e,om,Dth)
%
%       Disegna l'orbita nel piano perifocale
%
%       Dth è un vettore 1x2 contenente la theta iniziale e finale
%       dell'orbita da disegnare; funziona solo sulle ellissi
%
%

if(nargin==3)
    Dth=[0,2*pi];   % Se Dth non immesso fa il plot di tutta l'ellise
end

prec=150; % Indica la precisione (num di punti interpolati)

p=a*(1-e^2);

if(e>=1)
    c=a*e;
    b=sqrt(c^2-a^2);
    alpha=atan(b/a);
    theta=linspace(alpha-pi+0.1,pi-alpha-0.1,prec);
else
    theta=linspace(Dth(1),Dth(2),prec);
end

Rx=zeros(1,prec);
Ry=zeros(1,prec);
i=1;
for th=theta(1:prec)
    R=p/(1+e*cos(th));
    Rx(i)=R*cos(th+om);
    Ry(i)=R*sin(th+om);
    i=i+1;
end

plot(Rx,Ry,'LineWidth',1.5);
hold on;
axis equal;
