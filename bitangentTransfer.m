function [deltaV1, deltaV2, deltaV, deltat] = bitangentTransfer(ai, ei, af, ef, type, mu)
%{
INPUT: 
        - orbIniz:  [6x1], parametri orbita iniziale
        - orbFin:  [6x1], parametri orbita iniziali, theta deve rimanere
        uguale
        - type: sceglie il caso: 'ap', 'pa', 'pp', 'aa'
OUTPUT
        -Delta V1
        -Delta V2
        -Delta Vtot
        -Tempo di manovra (Delta t)
%}

rAi = ai*(1+ei);        %Calcolo raggi di apocentro e pericentro delle orbite fornite
rPi = ai*(1-ei);
rAf = af*(1+ef);
rPf = af*(1-ef);

switch lower(type)
    case 'ap'                           %Trasferimento apocentro orbita 1 - pericentro orbita 2
        at = (rAi+rPf)/2;
        et = abs((rAi-rPf)/(rAi+rPf));
        
        deltaV1 = abs(sqrt(2*mu*((1/rAi)-(1/(2*at)))) - sqrt(2*mu*((1/rAi)-(1/(2*ai)))));
        deltaV2 = abs(sqrt(2*mu*((1/rPf)-(1/(2*af)))) - sqrt(2*mu*((1/rPf)-(1/(2*at)))));
        
        deltaV = deltaV1 + deltaV2;
        
        deltat = pi() * sqrt((at^3)/mu);
        
    case 'pa'                           %Trasferimento pericentro orbita 1 - apocentro orbita 2
        at = (rAf+rPi)/2;
        et = abs((rAf-rPi)/(rAf+rPi));
        
        deltaV1 = abs(sqrt(2*mu*((1/rPi)-(1/(2*at)))) - sqrt(2*mu*((1/rPi)-(1/(2*ai)))));
        deltaV2 = abs(sqrt(2*mu*((1/rAf)-(1/(2*af)))) - sqrt(2*mu*((1/rAf)-(1/(2*at)))));
        
        deltaV = deltaV1 + deltaV2;
        
        deltat = pi() * sqrt((at^3)/mu);
    
    case 'aa'                           % Trasferimento da un apocentro all'altro
        at = (rPi + rPf)/2;
        et = abs((rAf-rAi)/(rAf + rAi));
        
        deltaV1 = abs(sqrt(2*mu*((1/rAi)-(1/(2*at)))) - sqrt(2*mu*((1/rAi)-(1/(2*ai)))));
        deltaV2 = abs(sqrt(2*mu*((1/rAf)-(1/(2*af)))) - sqrt(2*mu*((1/rAf)-(1/(2*at)))));
        
        deltaV = deltaV1 + deltaV2;
        
        deltat = pi() * sqrt((at^3)/mu);
        
    case 'pp'                           %Trasferimento da un pericentro all'altro
        
        at = (rPf + rPi)/2;
        et = abs((rPf - rPa)/(rPf + rPa));
        
        deltaV1 = abs(sqrt(2*mu*((1/rPi)-(1/(2*at)))) - sqrt(2*mu*((1/rPi)-(1/(2*ai)))));
        deltaV2 = abs(sqrt(2*mu*((1/rPf)-(1/(2*af)))) - sqrt(2*mu*((1/rPf)-(1/(2*at)))));
        
        deltaV = deltaV1 + deltaV2;
        
        deltat = pi() * sqrt((at^3)/mu);
end

