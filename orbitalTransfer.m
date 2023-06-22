function [DvA,DvB,th,Dt,P] = orbitalTransfer (type,a,e,om,rt)
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
%           P = parametri orbita/orbite di trasferimento; se le orbite di
%               trasferimento sono più di una la seconda viene messa nella
%               seconda riga.
%               Vengono date solo le orbite del caso più vantaggioso in
%               termini di costo. **Tranquillamente modificabile
%               Es: P = [a_t e_t]
%               Per 'circ' è il raggio dell'orbita circolarizzata
%               Es : P = [at_aa(1) et_aa(1) at_aa(2) et_aa(2);
%                         at_pp(1) et_pp(1) at_pp(2) et_pp(2)];
%                   Per biell
%               *******IMPLEMENTATO SOLO IN bitan e circ e biell e hohmann PER ORA*******
%
%-------------------------------------TYPE----------------------------------%
%
%           rotom = rotazione dell'orientazione dell'orbita;
%                   mantiene l'orbita invariata tranne per il parametro om
%           circ = manovra che circolarizza in un punto theta assegnato
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
%               a = semiasse maggiore iniziale e finale
%               1x1 per il caso 'rotom' e 'circ'
%               1x2 per gli alti casi
%
%           e : Vettore di dimensione 1x1 o 1x2; 
%               1x1 per il caso 'rotom', 'hohmann' e 'circ':
%               e = eccentricità iniziale e finale
%               1x2 per gli altri casi
%                   e(1) = eccentricità iniziale
%                   e(2) = eccentricità finale
%
%           om : Vettore di dimensione 1x2; [in radianti]
%               1x2 per il caso 'rotom', 'bitan', 'biell' e 'istant':
%                   om(1) = omega iniziale
%                   om(2) = omega finale
%               Per gli altri casi si può impostare anche nullo
%               NOTA: per 'biell' e 'bitan' basta solo che om siano uguali
%                     o shiftati di pi
%
%           rt : Caso 'biell' : raggio apocentro trasferimento biellittico;
%                Caso 'circ' : theta al quale far avvenire la
%                              circolarizzazione
%

mu=398600;

switch type
    case 'rotom'
        % Restituisce theta iniziale e finale sia di a che di b

        % Assegnazione di un singolo valore se vettore simmetrico:
        if size(a,2)==2
            if(a(1)==a(2))
                a=a(1);
            end
        end
        if size(e,2)==2
            if(e(1)==e(2))
                e=e(1);
            end
        end

        if(size(a,2)~=1 || size(e,2)~=1)
            fprintf("Errore inserimento dati (e o a)\n")
        end

        p=a*(1-e^2);    % p dell'orbita

        Dom=om(2)-om(1);

        % Se i pericentri puntano in direzioni opposte conviene girare il
        % minimo e quindi girare fino ad avere pericentri allineati con
        % versi opposti:
        
        if(Dom<-pi/2 || Dom>pi/2)
            Dom=Dom+pi;             % Da verificare
            Dom=minp2p(Dom);
            fprintf("Vero omega finale: %d \n",o2pi(om(2)+pi));
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

    case 'circ'

        % Assegnamento valori
        th = rt;
        a_i =  a;
        e_i = e;
        p_i = (a_i*(1-(e_i^2)));
        a_f = p_i/ (1+e_i*cos(th));

        % Calcolo costo
        v_tan_i = sqrt(mu/p_i)*(1+e_i*cos(th));
        v_tan_f = sqrt(mu/a_f);
        Dv_tang = abs(v_tan_f - v_tan_i);
        Dv_rad = abs( sqrt(mu/p_i)*e_i*sin(th));

        DvA = sqrt(Dv_rad^2 + Dv_tang^2);
        P = a_f;

        % Tempo = istantaneo e solo un possibile costo
        Dt=0;
        DvB=0;

    case 'hohmann'
        
        % Verifica dei valori immessi:
        if(size(a,2)~=2)
            fprintf("Errore nei valori immessi (dim a o e)\n")
        end

        v_i=sqrt(mu/a(1));      % Velocità nell'orbita iniziale
        v_f=sqrt(mu/a(2));      % Velocità nell'orbita finale

        a_t=(a(1)+a(2))/2;      % a dell'orbita di trasferimento
        e_t=abs( (a(1)-a(2)) / (a(1)+a(2)) );

        vp_t=sqrt(2*mu*((1/a(1))-(1/(2*a_t))));   % vel peric orbita trasf.
        va_t=sqrt(2*mu*((1/a(2))-(1/(2*a_t))));   % vel apoc orbita trasf.

        DvA=[abs(vp_t-v_i),abs(v_f-va_t)];   % Velocità di trasferimento

        T=2*pi*sqrt((a_t^3)/mu);
        Dt=T/2;                 % Tempo impiegato per la manovra

        DvB=0;                  % Solo per indicare che non esiste
        th=0;

        P = [a_t,e_t];

    case 'bitan' % Da controllare
        
        % Controllo valori inseriti
        if(size(a,2)~=2 || size(e,2)~=2)
            fprintf("Errore nell'inserimento dei valori\n")
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
                e_t_ap=abs((ra(1)-rp(2)) / (ra(1)+rp(2)));

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
                e_t_pa=abs((rp(1)-ra(2)) / (rp(1)+ra(2)));

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
                    P=[a_t_ap, e_t_ap;
                       a_t_pa, e_t_pa];
                    fprintf("Caso più conveniente: [A->P']\n");
                elseif Dvtot_ap>Dvtot_pa
                    DvA=[Dv1_pa,Dv2_pa];
                    DvB=[Dv1_ap,Dv2_ap];
                    Dt=[Dt_pa,Dt_ap];
                    th=[0,pi];
                    P=[a_t_pa, e_t_pa;
                       a_t_ap, e_t_ap];
                    fprintf("Caso più conveniente: [P->A']\n");
                end
                % NOTA: volendo puoi fare in modo che vengano restituiti
                %       tutti i valori utilizzando Dv1, ecc.. come vettori



            case pi

                %---------- Caso [A->A']: ----------%

                % Calcolo parametri orbita di trasferimento:
                a_t_aa=(ra(1)+ra(2))/2;
                e_t_aa=abs((ra(1)-ra(2)) / (ra(1)+ra(2)));

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
                e_t_pp=abs((rp(1)-rp(2)) / (rp(1)+rp(2)));

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
                    P=[a_t_aa, e_t_aa;
                       a_t_pp, e_t_pp];
                    fprintf("Caso più conveniente: [A->A']\n");
                elseif Dvtot_aa>Dvtot_pp
                    DvA=[Dv1_pp,Dv2_pp];
                    DvB=[Dv1_aa,Dv2_aa];
                    Dt=[Dt_pp,Dt_aa];
                    th=[0,0];
                    P=[a_t_pp, e_t_pp;
                       a_t_aa, e_t_aa];
                    fprintf("Caso più conveniente: [P->P']\n");
                end
                % NOTA: volendo puoi fare in modo che vengano restituiti
                %       tutti i valori utilizzando Dv1, ecc.. come vettori

            otherwise
                fprintf("Errore inserimento dati (le eccentricità devono avere la " + ...
                    "stessa direzione)\n");
        end


    case 'biell' % Da controllare

        % Controllo valori inseriti
        if(size(a,2)~=2 || size(e,2)~=2)
            fprintf("Errore nell'inserimento dei valori (dim a o e)\n")
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
                et_aa(1)=abs((ra(1)-rt)/(rt+ra(1)));
                et_aa(2)=abs((ra(2)-rt)/(rt+ra(2)));

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
                et_pp(1)=abs((rp(1)-rt)/(rt+rp(1)));
                et_pp(2)=abs((rp(2)-rt)/(rt+rp(2)));

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

                if  Dvtot_pp>Dvtot_aa
                    DvA=[Dv1_aa,Dv2_aa,Dv3_aa];
                    DvB=[Dv1_pp,Dv2_pp,Dv3_pp];
                    Dt=[Dt_aa,Dt_pp];
                    th=[pi,pi];
                    P = [at_aa(1) et_aa(1) at_aa(2) et_aa(2);
                         at_pp(1) et_pp(1) at_pp(2) et_pp(2)];
                    fprintf("Caso più conveniente: [A->A']\n");
                elseif Dvtot_aa>Dvtot_pp
                    DvA=[Dv1_pp,Dv2_pp,Dv3_pp];
                    DvB=[Dv1_aa,Dv2_aa,Dv3_aa];
                    Dt=[Dt_pp,Dt_aa];
                    th=[0,0];
                    P = [at_pp(1) et_pp(1) at_pp(2) et_pp(2);
                         at_aa(1) et_aa(1) at_aa(2) et_aa(2)];
                    fprintf("Caso più conveniente: [P->P']\n");
                end


            case pi

                %---------- Caso [A->P']: ----------%

                % Calcolo parametri orbite di trasferimento
                at_ap(1)=(ra(1)+rt)/2;
                at_ap(2)=(rt+rp(2))/2;
                et_ap(1)=abs((ra(1)-rt)/(rt+ra(1)));
                et_ap(2)=abs((rp(2)-rt)/(rt+rp(2)));

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
                et_pa(1)=abs((rp(1)-rt)/(rt+rp(1)));
                et_pa(2)=abs((ra(2)-rt)/(rt+ra(2)));

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
                    P = [at_ap(1) et_ap(1) at_ap(2) et_ap(2);
                         at_pa(1) et_pa(1) at_pa(2) et_pa(2)];
                    fprintf("Caso più conveniente: [A->P']\n");
                elseif Dvtot_ap>Dvtot_pa
                    DvA=[Dv1_pa,Dv2_pa,Dv3_pa];
                    DvB=[Dv1_ap,Dv2_ap,Dv3_ap];
                    Dt=[Dt_pa,Dt_ap];
                    th=[0,pi];
                    P = [at_pa(1) et_pa(1) at_pa(2) et_pa(2);
                         at_ap(1) et_ap(1) at_ap(2) et_ap(2)];
                    fprintf("Caso più conveniente: [P->A']\n");
                end


            
            otherwise
                fprintf("Errore inserimento dati (le eccentricità devono avere la " + ...
                    "stessa direzione)\n");
        end


        case 'istant' %%%%%%%% DA CONTROLLARE %%%%%%%%
            %%%%%%%%%%% ANZI DA RIFARE UTILIZZANDO POSZIONE E VELOCITA' AL
            %%%%%%%%%%% POSTO DEI PARAMETRI ORBITALI %%%%%%%%%%
            % Controllo dei valori immessi:
            if(size(a,2)~=2 || size(e,2)~=2)
                    fprintf("Errore nell'inserimento dei valori (dim a o e)\n")
            end


            %------------CALCOLO PUNTI DI INTERSEZIONE-----------%

            Dom = om(2)-om(1);
            prec = 100000;
            theta = linspace(0,2*pi,prec);
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

            %------------CALCOLO DEI COSTI PER OGNI PUNTO DI INTERSEZIONE-----------%

            % costo nel primo punto:
            [~, VV1] = paraorb2rv(a_t1,e_t1,i_t1,OM_t1,om_t1,inters(1)+Dom);
            [~, VV2] = paraorb2rv(a_t2,e_t2,i_t2,OM_t2,om_t2,inters(1));
            Dv_a = norm(VV1-VV2);

            % costo nel secondo punto:
            [~, VV1] = paraorb2rv(a_t1,e_t1,i_t1,OM_t1,om_t1,inters(2)+Dom);
            [~, VV2] = paraorb2rv(a_t2,e_t2,i_t2,OM_t2,om_t2,inters(2));
            Dv_b = norm(VV1-VV2);

            %%%%%% DA FINIRE: bisogna fare la scelta di quale punto usare e
            %%%%%% inserire i valori di output

    otherwise
        fprintf("Errore nella scelta del tipo\n");
end





