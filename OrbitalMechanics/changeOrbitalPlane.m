function [Dv,om_f,th]= changeOrbitalPlane (a,e,i_i,OM_i,om_i,i_f,OM_f,mu)
%
%       th indica la posizione (sia iniziale che finale) del punto di
%       manovra sia nell'orbita iniziale che quella finale
%
%       a,e,th Parametri fissi
%       i,OM,om Paramatri che definiscono l'orbita iniziale
%       Dv Parametro di cambio orbita
%       
%       Tutti i parametri angolari vanno inseriti e vengono restituiti in
%       radianti
%       a in km
%

if(nargin==7)
    mu=398600;
end

DOM=OM_f-OM_i;          % Delta omega delle due orbite
Di=i_f-i_i;             % Delta inclinazione delle due orbite

alpha=cos(i_i)*cos(i_f)+sin(i_i)*sin(i_f)*cos(DOM);
alpha=acos(alpha);

%{
if(sign(alpha)~=sign(DOM))
    alpha=-alpha;       % Inverte segno alpha se necessario
end
%}

if (alpha==0)
    fprintf('Stessa orbita \n');
end


% Calcolo per ogni caso il coseno di u_i e u_f e th e om_f:

if(DOM>0 && Di>0)
    Cu_i=(cos(i_f)-cos(alpha)*cos(i_i))/(-sin(alpha)*sin(i_i));
    Cu_f=(cos(i_i)-cos(alpha)*cos(i_f))/(sin(alpha)*sin(i_f));
end
if(DOM<0 && Di<0)
    Cu_i=(cos(i_f)-cos(alpha)*cos(i_i))/(-sin(alpha)*sin(i_i));
    Cu_f=(cos(i_i)-cos(alpha)*cos(i_f))/(sin(alpha)*sin(i_f));
end
if(DOM>0 && Di<0)
    Cu_i=(cos(i_f)-cos(alpha)*cos(i_i))/(sin(alpha)*sin(i_i));
    Cu_f=(cos(i_i)-cos(alpha)*cos(i_f))/-(sin(alpha)*sin(i_f));
end
if(DOM<0 && Di>0)
    Cu_i=(cos(i_f)-cos(alpha)*cos(i_i))/(sin(alpha)*sin(i_i));
    Cu_f=(cos(i_i)-cos(alpha)*cos(i_f))/-(sin(alpha)*sin(i_f));
end
if(DOM==0)
    Cu_i=1;
    Cu_f=1;
end


% Calcolo seno di i_i e i_f
Su_i=(sin(DOM)/sin(alpha))*sin(i_f);
Su_f=(sin(DOM)/sin(alpha))*sin(i_i);


% Calcolo u_i e u_f come tan(sen/cos) ma utilizzando atan2 che li riporta
% tra -pi e +pi
u_i=atan2(Su_i,Cu_i);
u_f=atan2(Su_f,Cu_f);

% Print utile per il controllo:
%fprintf('\n sen(u_i)=%f \n cos(u_i)=%f \n sen(u_f)=%f \n cos(u_f)=%f \n u_i=%f \n u_f=%f\n', Su_i, Cu_i, Su_f, Cu_f, u_i, u_f);


% Se delta omega nullo:
if(DOM==0)
    u_i=0;
    u_f=0;
end



if(DOM*Di>0 || (DOM>0 && Di==0))
    th=u_i-om_i;
    th=o2pi(th);
    om_f=u_f-th;
end
if(DOM*Di<0 || (DOM<0 && Di==0))
    th=2*pi-om_i-u_i;
    th=o2pi(th);
    om_f=2*pi-th-u_f;
end
if(DOM==0)
    th=-om_i;
    om_f=om_i;
end

om_f=o2pi(om_f);

if(th<pi/2 || th>(3*pi)/2)
    th = th+pi;         % nel caso theta sia nel primo quarto quadrante
end                     % conviene manovrare nell'altro punto di
                        % intersezione

% Nel caso si puo' verificare che v_theta con questa theta sia minore
% Oppure si puo' verificare con il segno del coseno 

v_th=sqrt(mu/(a*(1-e^2)))*(1+e*cos(th));    % Calcolo v tangente in P

Dv=2*v_th*sin(alpha/2);                     % Calcolo costo manovra 




