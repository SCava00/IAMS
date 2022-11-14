function [Dt] = deltaTime(a,e,th)
%
%       Funzione che ricevuto in ingresso i parametri di un orbita,
%       e le cordinate di due o un punto, restituisce il tempo di
%       percorrenza da l primo punto al secondo
%
%       Se th è uno scalare restituisce il tempo impiegato rispetto al
%       pericentro
%       Es: th = [ th_i , th_f ]
%       NOTA: -th va espresso in radianti
%             -a deve essere negativo per orbite iperboliche
%

% Verifica del tipo di orbita:
if(e<1)
    par=1;      % Orbite ellittiche
end
if(e==1)
    par=2;      % Orbite paraboliche
end
if(e>1)
    par=3;      % Orbite iperboliche
end


% Se th scalare imposta punto iniziale al pericentro
if(size(th,2)==1)
    th=[0,th];
end


% Funzioni che trasportano i valori di th tra 0 e 2pi radianti
th(1)=o2pi(th(1));
th(2)=o2pi(th(2));




switch par

    case 1          % Orbite ellittiche
        t=[0,0];
        for i=1:2
            if (th(i)==pi || th(i)==0)      
                t(i)=th(i)*sqrt(( a^3) / mu);
            else
                E=2*atan(tan(th(i)/2) * sqrt( (1-e)/(1+e) ) );
                t(i)=sqrt( (a^3)/mu ) * (E-e*sin(E));
            end
        end

        if(t(2)>=t(1))
            Dt=t(2)-t(1);
        end
        if(t(1)>t(2))
            T=2*pi*sqrt((a^3)/mu);
            Dt=T-(t(1)-t(2));
        end

    case 2          % Orbite paraboliche

        if(a>0)
            a=-a;
        end

        p=-a*(1-e^2);

        % Calcolo dell'angolo limite a cui può arrivare
        c=a*e;
        b=sqrt(c^2-a^2);
        alpha=-atan(b/a);

        D=[0,0];
        t=[0,0];

        for i=1:2
            if(th(i)>pi-alpha || th(i)<pi+alpha)
                fprintf("Theta non esistente nell'orbita");
            else
                D(i)=tan(th(i)/2);
                t(i)=0.5*sqrt((p)^3/mu)*(D(i) + ((D(i)^3)/3) );
            end
        end

        Dt=abs(t(2)-t(1));



    case 3          % Orbite iperboliche

        if(a>0)
            a=-a;
        end

        % Calcolo dell'angolo limite a cui può arrivare
        c=a*e;
        b=sqrt(c^2-a^2);
        alpha=-atan(b/a);

        F=[0,0];
        t=[0,0];

        for i=1:2
            if(th(i)>pi-alpha || th(i)<pi+alpha)
                fprintf("Theta non esistente nell'orbita");
            else
                F(i)=2*atanh( sqrt( (e-1)/(e+1) ) * tan(th(i)/2) );
                t(i)=sqrt(((-a)^3)/mu)*( e*sinh(F) - F);
            end
        end

        Dt=abs(t(2)-t(1));
end

