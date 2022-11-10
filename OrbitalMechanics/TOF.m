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


E1 = 2 * atan(sqrt((1 - e)/(1 + e)) * tan(th1/2));                     %Calcolo anomalia eccentrica E1
E2 = 2 * atan(sqrt((1 - e)/(1 + e)) * tan(th2/2));                     %Calcolo anomalia eccentrica E2

deltat = sqrt((a^3)/mu) * ((E2-E1) - e * (sin(E2) - sin(E1)));         %Calcolo TOF

if th2<th1                                                             %
    deltat = deltat + 2 * pi() * sqrt((a^3) / mu);

end

