function [deltat] = TOF(a, e, th1, th2, mu)
%{  
    Funzione Time of Flight

    INPUT 
         - Semiasse maggiore dell'orbita
         - Eccentricità dell'orbita
         - Anomalia vera iniziale
         - Anomalia vera finale
         - Costante gravitazionale planetaria mu in [km3/s2] (per la terra inserire
         398600)
    OUTPUT
         - Tempo di trasferimento tra le anomalie
%}

dth = th2-th1;

E = 2 * atan(sqrt((1 - e)/(1 + e)) * tan(dth/2)); %Calcolo anomalia eccentrica E

M = E - e*sin(E);                                 %Calcolo anomalia media M

deltat = sqrt((a^3)/mu) * M;                      %Calcolo TOF

end

