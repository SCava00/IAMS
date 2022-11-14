function [DvA,DvB,th,Dt] = orbitalTransfer (type,a,e,om,rt)
%
%
%------------------------------------OUTPUT---------------------------------%
%
%           DvA = Vettore 1x1, 1x2 o 1x3 in base al tipo selezionato.
%                 Contiene le variazioni di velocità del cambio orbitale
%                 più conveniente (in termini di costo totale).
%                 Es: DvA = [Dv1,Dv2,Dv3]
%
%           DvB = Vettore 1x2 o 1x3 in base al tipo selezionato.
%                 Contiene le variazioni di velocità del cambio orbitale 
%                 meno conveniente (in termini di costo totale).
%                 NOTA : Vettore non nullo solo nel caso 'biell' e 'bitan'
%
%           th = Caso 'rotom':
%                     th è una matrice 2x2 che contine i theta iniziali e
%                     finali dei due punti in cui si può fare avvenire la
%                     manovra
%                         |tha_i    tha_f|
%                     th= |              |
%                         |thb_i    thb_f|
%
%                Caso 'bitan' e 'biell':
%                     th è un vettore 1x2 che indica la theta dove avviene
%                     la prima manovra sulla orbita iniziale e la theta dove
%                     avviene l'ultima manovra, sulla orbita finale
%                     La manovra considerata è quella più conveniente in
%                     termini di costo totale; è qindi quella di DvA
%                     In pratica 0 = pericentro
%                                pi = apocentro
%                     Es: th = [0,pi] significa che la manovra più
%                         conveniente è quella del tipo [P->A']
%
%                Caso 'hohmann':
%                     th è nullo perchè inutile
%
%                Caso 'istant':
%                     th è una matrice 2x2 che contiene i valori di theta
%                     iniziali e finali dei due punti di intersezione delle
%                     due orbite:
%                         |tha_i    tha_f|
%                     th= |              |
%                         |thb_i    thb_f|
%
%                 
%           Dt = Vettore 1x1 o 1x2 a seconda del tipo scelto.
%                Nel caso 'rotom' e 'istant'  è nullo perchè istantaneo.
%                Nel caso 'hohmann' è il tempo impiegato per cambiare
%                orbita.
%                Nei casi 'biell' e 'bitan' è un vettore contenete il tempo
%                impegato nel trasferimento orbitale più economico e in
%                quello meno economico in termini di costo totale.
%                Es: Dt = [ Dt(DvA), Dt(DvB) ]
%
%
%
%-------------------------------------TYPE----------------------------------%
%
%           rotom = rotazione dell'orientazione dell'orbita;
%                   mantiene l'orbita invariata tranne per il parametro om
%           hohmann = manovra alla hohmann (orbite circolari); in questo
%                     caso Dv3=0 
%           bitan = trasferimento ellittico bitangente, condizione
%                    necessaria è che i vettori eccentricità abbiano stessa
%                    direzione; Dv3=0
%           biell = trasferimento biellittico
%           istant = trasferimento istantaneo nei punti di intersezone
%
%------------------------------------INPUT----------------------------------%
%
%           a : Vettore di dimensione 1x1 o 1x2:
%               1x1 per il caso 'rotom' e 'hohmann'
%               1x2 per gli alti casi
%
%           e : Vettore di dimensione 1x1 o 1x2; 
%               1x1 per il caso 'rotom' e 'hohmann':
%               e = eccentricità iniziale e finale
%               1x2 per gli altri casi
%                   e(1) = eccentricità iniziale
%                   e(2) = eccentricità finale
%
%           om : Vettore di dimensione 1x2; [in radianti]
%               1x2 per il caso 'rotom', 'bitan' e 'biell':
%                   om(1) = omega iniziale
%                   om(2) = omega finale
%               Per gli altri casi si può impostare anche nullo
%               NOTA: per 'biell' e 'bitan' basta solo che om siano uguali
%                     o shiftati di pi
%
%           rt : Raggio apocentro trasferimento biellittico;
%                utilizzato solo nel caso 'biell'
%

% Test:
%{
type='bitang';
a=[10000, 15000];
e=[0.1,0.1];
om=[0,0];
%}

mu=398600;

switch type
    case 'rotom'
        % Restituisce theta iniziale e finale sia di a che di b

        % Assegnazione di un singolo valore se vettore simmetrico:
        if(a(1)==a(2))
            a=a(1);
        end
        if(e(1)==e(2))
            e=e(1);
        end

        if(size(a,2)~=1 || size(e,2)~=1)
            fprintf("Errore inserimento dati (e o a)")
        end

        p=a*(1-e^2);    % p dell'orbita

        Dom=om(2)-om(1);

        % Se i pericentri puntano in direzioni opposte conviene girare il
        % minimo e quindi girare fino ad avere pericentri allineati con
        % versi opposti:
        
        if(Dom<-pi/2 || Dom>pi/2)
            Dom=Dom+pi;             % Da verificare
            Dom=minp2p(Dom);
            fprintf("Vero omega finale: %d",o2pi(om(2)+pi));
        end
        

        tha_i=Dom/2;              % Calcolo dei punti di manovra 
        thb_i=tha_i+pi;           % rispetto alla orbita iniziale

        tha_f=2*pi-tha_i;         % Calcolo punti di manovra 
        thb_f=2*pi-thb_i;         % rispetto orbita finale

        DvA=2*sqrt(mu/p)*e*sin(Dom/2);  % Calcolo Dv risultante

        th = [tha_i,tha_f;thb_i,thb_f];

        Dt=0;           % Nel senso che è istantaneo appena si raggiunge 
                        % il punto di manovra

        DvB=0;          % Solo per indicare che non esiste

    case 'hohmann'

        % Assegnazione di un singolo valore se vettore simmetrico:
        if(a(1)==a(2))
            a=a(1);
        end
        if(e(1)==e(2))
            e=e(1);
        end
        
        % Verifica dei valori immessi:
        if(e~=1 || size(e,2)~=1  || size(a,2)~=2)
            fprintf("Errore nei valori immessi (dim a o e)")
        end

        v_i=sqrt(mu/a(1));      % Velocità nell'orbita iniziale
        v_f=sqrt(mu/a(2));      % Velocità nell'orbita finale

        a_t=(a(1)+a(2))/2;      % a dell'orbita di trasferimento

        vp_t=sqrt(2*mu*((1/a(1))-(1/(2*a_t))));   % vel peric orbita trasf.
        va_t=sqrt(2*mu*((1/a(2))-(1/(2*a_t))));   % vel apoc orbita trasf.

        DvA=[abs(vp_t-v_i),abs(v_f-va_t)];   % Velocità di trasferimento

        T=2*pi*sqrt((a_t^3)/mu);
        Dt=T/2;                 % Tempo impiegato per la manovra

        DvB=0;                  % Solo per indicare che non esiste
        th=0;


    case 'bitan' % Da controllare
        
        % Controllo valori inseriti
        if(size(a,2)~=2 || size(e,2)~=2)
            fprintf("Errore nell'inserimento dei valori")
        end


        % Calcolo raggi apocentro e pericentro delle due orbite
        ra(1)=a(1)*(1+e(1));
        ra(2)=a(2)*(1+e(2));
        rp(1)=2*a(1)-ra(1);
        rp(2)=2*a(2)-ra(2);


        Dom=abs(om(1)-om(2));        % Calcolo angolo tra i vettori eccentricità

        switch Dom
            case 0

                %---------- Caso [A->P']: ----------%

                % Calcolo parametri orbita di trasferimento:
                a_t_ap=(ra(1)+rp(2))/2;

                % Calcolo Dv1:
                Dv1_ap=deltaVtang (a(1),a_t_ap,ra(1));

                % Calcolo Dv2:
                Dv2_ap=deltaVtang (a_t_ap,a(2),rp(2));

                % Calcolo Dvtot caso [A->P']:
                Dvtot_ap=Dv1_ap+Dv2_ap;

                % Calcolo Dt caso [A->P']:
                Dt_ap=pi*sqrt((a_t_ap^3)/mu);

             

                %---------- Caso [P->A']: ----------%

                % Calcolo parametri orbita di trasferimento:
                a_t_pa=(rp(1)+ra(2))/2;

                % Calcolo Dv1:
                Dv1_pa= deltaVtang (a(1),a_t_pa,rp(1));

                % Calcolo Dv2
                Dv2_pa= deltaVtang (a_t_pa,a(2),ra(2));

                % Calcolo Dvtot caso [P->A']:
                Dvtot_pa=Dv1_pa+Dv2_pa;

                % Calcolo Dt caso [P->A']:
                Dt_pa=pi*sqrt((a_t_pa^3)/mu);

                

                % Verifica di quale sia più conveniente:

                if Dvtot_pa>Dvtot_ap
                    DvA=[Dv1_ap,Dv2_ap];
                    DvB=[Dv1_pa,Dv2_pa];
                    Dt=[Dt_ap,Dt_pa];
                    th=[pi,0];
                    fprintf("Caso più conveniente: [A->P']");
                elseif Dvtot_ap>Dvtot_pa
                    DvA=[Dv1_pa,Dv2_pa];
                    DvB=[Dv1_ap,Dv2_ap];
                    Dt=[Dt_pa,Dt_ap];
                    th=[pi,0];
                    fprintf("Caso più conveniente: [P->A']");
                end
                % NOTA: volendo puoi fare in modo che vengano restituiti
                %       tutti i valori utilizzando Dv1, ecc.. come vettori



            case pi

                %---------- Caso [A->A']: ----------%

                % Calcolo parametri orbita di trasferimento:
                a_t_aa=(ra(1)+ra(2))/2;

                % Calcolo Dv1:
                Dv1_aa=deltaVtang (a(1),a_t_aa,ra(1));

                % Calcolo Dv2:
                Dv2_aa= deltaVtang (a_t_aa,a(2),ra(2));

                % Calcolo Dvtot caso [A->A']:
                Dvtot_aa=Dv1_aa+Dv2_aa;

                % Calcolo Dt caso [A->A']:
                Dt_aa=pi*sqrt((a_t_aa^3)/mu);

             

                %---------- Caso [P->P']: ----------%

                % Calcolo parametri orbita di trasferimento:
                a_t_pp=(rp(1)+rp(2))/2;

                % Calcolo Dv1:
                Dv1_pp= deltaVtang (a(1),a_t_pp,rp(1));

                % Calcolo Dv2
                Dv2_pp= deltaVtang (a_t_pp,a(2),rp(2));

                % Calcolo Dvtot caso [P->P']:
                Dvtot_pp=Dv1_pp+Dv2_pp;

                % Calcolo Dt caso [P->P']:
                Dt_pp=pi*sqrt((a_t_pp^3)/mu);

               

                % Verifica di quale sia più conveniente:

                if Dvtot_pp>Dvtot_aa
                    DvA=[Dv1_aa,Dv2_aa];
                    DvB=[Dv1_pp,Dv2_pp];
                    Dt=[Dt_aa,Dt_pp];
                    th=[pi,pi];
                    fprintf("Caso più conveniente: [A->A']");
                elseif Dvtot_aa>Dvtot_pp
                    DvA=[Dv1_pp,Dv2_pp];
                    DvB=[Dv1_aa,Dv2_aa];
                    Dt=[Dt_pp,Dt_aa];
                    th=[0,0];
                    fprintf("Caso più conveniente: [P->P']");
                end
                % NOTA: volendo puoi fare in modo che vengano restituiti
                %       tutti i valori utilizzando Dv1, ecc.. come vettori

            otherwise
                fprintf("Errore inserimento dati (le eccentricità devono avere la " + ...
                    "stessa direzione)");
        end


    case 'biell' % Da controllare

        % Controllo valori inseriti
        if(size(a,2)~=2 || size(e,2)~=2)
            fprintf("Errore nell'inserimento dei valori (dim a o e)")
        end


        % Calcolo raggi apocentro e pericentro delle due orbite
        ra(1)=a(1)*(1+e(1));
        ra(2)=a(2)*(1+e(2));
        rp(1)=2*a(1)-ra(1);
        rp(2)=2*a(2)-ra(2);


        Dom=abs(om(1)-om(2));        % Calcolo angolo tra i vettori eccentricità

        switch Dom
            case 0
                %---------- Caso [A->A']: ----------%

                % Calcolo parametri orbite di trasferimento
                at_aa(1)=(ra(1)+rt)/2;
                at_aa(2)=(rt+ra(2))/2;

                % Calcolo Dv1:
                Dv1_aa= deltaVtang (a(1),at_aa(1),ra(1));

                % Calcolo Dv2:
                Dv2_aa= deltaVtang (at_aa(1),at_aa(2),rt);

                % Calcolo Dv3:
                Dv3_aa= deltaVtang (at_aa(2),a(2),ra(2));

                % Calcolo Dt caso [A->A']:
                Dt_aa=pi*( sqrt((at_aa(1))/mu) + sqrt((at_aa(2))/mu));


                %---------- Caso [P->P']: ----------%

                % Calcolo parametri orbite di trasferimento:
                at_pp(1)=(rp(1)+rt)/2;
                at_pp(2)=(rt+rp(2))/2;

                % Calcolo Dv1
                Dv1_pp= deltaVtang (a(1),at_pp(1),rp(1));

                % Calcolo Dv2
                Dv2_pp= deltaVtang (at_pp(1),at_pp(2),rt);

                % Calcolo Dv3
                Dv3_pp= deltaVtang (at_pp(2),a(2),rp(2));

                % Calcolo Dt caso [P->P']:
                Dt_pp=pi*( sqrt((at_pp(1))/mu) + sqrt((at_pp(2))/mu));


                % Verifica di quale caso sia più conveniente:

                Dvtot_pp=Dv1_pp+Dv2_pp+Dv3_pp;
                Dvtot_aa=Dv1_aa+Dv2_aa+Dv3_aa;

                if Dvtot_pp>Dvtot_aa
                    DvA=[Dv1_aa,Dv2_aa,Dv3_aa];
                    DvB=[Dv1_pp,Dv2_pp,Dv3_pp];
                    Dt=[Dt_aa,Dt_pp];
                    th=[pi,pi];
                    fprintf("Caso più conveniente: [A->A']");
                elseif Dvtot_aa>Dvtot_pp
                    DvA=[Dv1_pp,Dv2_pp,Dv3_pp];
                    DvB=[Dv1_aa,Dv2_aa,Dv3_aa];
                    Dt=[Dt_pp,Dt_aa];
                    th=[0,0];
                    fprintf("Caso più conveniente: [P->P']");
                end


            case pi

                %---------- Caso [A->P']: ----------%

                % Calcolo parametri orbite di trasferimento
                at_ap(1)=(ra(1)+rt)/2;
                at_ap(2)=(rt+rp(2))/2;

                % Calcolo Dv1:
                Dv1_ap=deltaVtang (a(1),at_ap(1),ra(1));

                % Calcolo Dv2:
                Dv2_ap= deltaVtang (at_ap(1),at_ap(2),rt);

                % Calcolo Dv3:
                Dv3_ap= deltaVtang (at_ap(2),a(2),rp(2));

                % Calcolo Dt caso [A->P']:
                Dt_ap=pi*( sqrt((at_ap(1))/mu) + sqrt((at_ap(2))/mu));



                %---------- Caso [P->A']: ----------%

                % Calcolo parametri orbite di trasferimento:
                at_pa(1)=(rp(1)+rt)/2;
                at_pa(2)=(rt+ra(2))/2;

                % Calcolo Dv1
                Dv1_pa= deltaVtang (a(1),at_pa(1),rp(1));

                % Calcolo Dv2
                Dv2_pa= deltaVtang (at_pa(1),at_pa(2),rt);

                % Calcolo Dv3
                Dv3_pa= deltaVtang (at_pa(2),a(2),ra(2));

                % Calcolo Dt caso [P->A']:
                Dt_pa=pi*( sqrt((at_pa(1))/mu) + sqrt((at_pa(2))/mu));



                % Verifica di quale caso sia più conveniente:

                Dvtot_pa=Dv1_pa+Dv2_pa+Dv3_pa;
                Dvtot_ap=Dv1_ap+Dv2_ap+Dv3_ap;

                if Dvtot_pa>Dvtot_ap
                    DvA=[Dv1_ap,Dv2_ap,Dv3_ap];
                    DvB=[Dv1_pa,Dv2_pa,Dv3_pa];
                    Dt=[Dt_ap,Dt_pa];
                    th=[pi,0];
                    fprintf("Caso più conveniente: [A->P']");
                elseif Dvtot_ap>Dvtot_pa
                    DvA=[Dv1_pa,Dv2_pa,Dv3_pa];
                    DvB=[Dv1_ap,Dv2_ap,Dv3_ap];
                    Dt=[Dt_pa,Dt_ap];
                    th=[0,pi];
                    fprintf("Caso più conveniente: [P->A']");
                end


            
            otherwise
                fprintf("Errore inserimento dati (le eccentricità devono avere la " + ...
                    "stessa direzione)");
        end


        case 'istant'

            % Controllo dei valori immessi:
            if(size(a,2)~=2 || size(e,2)~=2)
                    fprintf("Errore nell'inserimento dei valori (dim a o e)")
            end

            if (a(1)<0)
                a(1)=-a(1);
            end
            if (a(2)>0)
                a(2)=-a(2);
            end

            %------------CALCOLO PUNTI DI INTERSEZIONE-----------%

            % Calcolo parametri orbite: 
            p(1)=a(1)*(1-e(1)^2);
            p(2)=a(2)*(1-e(2)^2);

            % Calcolo parametri di risoluzioni dei punti di intersezione
            Dom= om(2)-om(1);
            parA = p(2)*e(1)*cos(Dom) - p(1)*e(2);
            parB = -p(2)*e(1)*sin(Dom);
            parC = p(2)-p(1);

            % Definizione funzione non lineare da risolvere
            f = @(x) cos(x)*parA + sin(x)*parB + parC;

            % Risoluzione funzione (calcolo punti intersezione rispetto
            % alla seconda orbita)
            thA(2)=fzero(f,0);
            thB(2)=fzero(f,pi);

            % Calcolo punti di intersezione rispetto prima orbita
            thA(1)=thA(2)+pi;
            thB(2)=thB(2)+pi;



            %-------CALCOLO VELOCITA' NEI PUNTI DI INTERSEZIONE-------%

            par(1)=sqrt(mu/p(1)); % Ottimizzazione banale
            par(2)=sqrt(mu/p(2));

            % Velocità lungo la prima orbita
            vthA(1)=par(1)*(1+e(1)*cos(thA(1)));
            vthB(1)=par(1)*(1+e(1)*cos(thB(1)));

            vrA(1)=par(1)*e(1)*sin(thA(1));
            vrB(1)=par(1)*e(1)*sin(thB(1));


            % Velocità lungo la seconda orbita
            vthA(2)=par(2)*(1+e(2)*cos(thA(2)));
            vthB(2)=par(2)*(1+e(2)*cos(thB(2)));

            vrA(2)=par(2)*e(2)*sin(thA(2));
            vrB(2)=par(2)*e(2)*sin(thB(2));



            %----------CALCOLO DEI DELTA V NECESSARI----------%

            DvthA=vthA(2)-vthA(1);
            DvrA=vrA(2)-vrA(1);

            DvA=sqrt(DvthA^2 + DvrA);

            DvthB=vthB(2)-vthB(1);
            DvrB=vrB(2)-vrB(1);

            DvB=sqrt(DvthB^2 + DvrB^2);

            % Controllo di quale sia effettivamente maggiore, e
            % assegnamento dei valori di theta:
            if DvA>DvB
                par=DvA;
                DvA=DvB;
                DvB=par;

                th=[thB(1),thB(2);thA(1),thA(2)];
            else
                th=[thA(1),thA(2);thB(1),thB(2)];
            end

            Dt=0;           % Perhè istantaneo una volta arrivati a quel valore di theta



    otherwise
        fprintf("Errore nella scelta del tipo\n");
end





