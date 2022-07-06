%matriz H diferentes distancias

%Antenas receptoras (array de NR elementos)
NR = 16; %nomes NR/2 parells +1
f = 30e9; %antRx dipolo no permite f superiores, dificultad calculo
frx = 30e9;
c = 299792458;
lambda = c/f;
kl = 2*pi/lambda;
%D = 645; % distancia franhofer
%dRx = 21*lambda; %2lambda
y = 0;
VecDis = 100*1.259.^[0:1:13];

NF = 64; %64
NC = 64; %64
N = NF*NC;
d = 2*lambda; %2lambda
dA = 0;
dRx = 12*lambda;
Agrupacio = 16;
canvi_sector = 0;
dAnombre = 0;

D = 1e4;

AntenasRx = zeros(NR,3);
%AntenasRx(:,2) = y;
AntenasRx(:,3) = D; %(z)
inc_fase = -pi/4;
canvi_DQ = -0.25;

%Generació RIS N elements
RIS = zeros(NF, NC); %x
RIS (:,:,2) = zeros(NF, NC); %y
RIS (:,:,3) = zeros(NF, NC); %z
RIS(1,:,1) = -NC/2*d - NF/(2*Agrupacio)*dA +d*0.5+dA*0.5;

    for j = 2:NC
        canvi_sector = canvi_sector+1;
        if canvi_sector >= Agrupacio 
            RIS(j,:,1) = d + RIS(j-1,1,1) +dA;
            canvi_sector = 0;
        else 
            RIS(j,:,1) = d + RIS(j-1,1,1);
        end
    end
    
RIS(:,1,2) = -NF/2*d - NF/(2*Agrupacio)*dA +d*0.5+dA*0.5 ;

    for j = 2:NF
        canvi_sector = canvi_sector+1;
        if canvi_sector >= Agrupacio 
            RIS(:,j,2) = d + RIS(1,j-1,2) +dA;
            canvi_sector = 0;
        else 
            RIS(:,j,2) = d + RIS(1,j-1,2);
        end
    end

RISv = zeros(N,3);
RISv(:,1) = reshape(RIS(:,:,1), N,1);
RISv(:,2) = reshape(RIS(:,:,2), N,1);
RISv(:,3) = reshape(RIS(:,:,3), N,1);

%Generacio RX NR elements
DistQuadrados = RIS(1,1,1)/(NR/4);
%AntenasRx(1,1) = 0;

AntenasRx(1,1) = -(0.5*dRx+3*dRx);
AntenasRx(NR/2+1,2) = -(0.5*dRx+3*dRx);

    for j = 2:NR/2
        AntenasRx(j,1) = AntenasRx(j-1,1)+dRx;
    end
    for j= NR/2+2:NR
        AntenasRx(j,2) = AntenasRx(j-1,2)+dRx;
    end

%Calcul matriu Distancies, Azimut i elevació de element i-essim respecte
%recetptor i-essim de AntenasRx. Calcul matriu Directivitat
Distancies = zeros(NR,N);
Azimut = zeros(NR,N);
Elevacio = zeros(NR,N);
AzimutRx = zeros(NR,N);
ElevacioRx = zeros(NR,N);

Directivitat = zeros(NR, N);
faseAntena = zeros(NR, N);eTx = zeros(N*NR,3); %depenent de la pol -> ,1 o ,3

DirectivitatRx = zeros(NR, N);
faseAntenaRx = zeros(NR, N);
eRx = zeros(N*NR,3); %depenent de la pol -> ,1 o ,3

Antenna = patchMicrostrip('Height' ,9.9931e-05,'Length',0.0048, 'Width',0.0062,'GroundPlaneLength', 0.0100, 'GroundPlaneWidth', 0.0100, 'PatchCenterOffset', [0 0],'FeedOffset', [0.0010 0],'Load', lumpedElement );
antRx = dipole('Length', 0.0046967 , 'Width', 9.9931e-05, 'Tilt',90,'TiltAxis',[1 0 0]);
%antRx = patchMicrostrip('Height' ,9.9931e-05,'Length',0.0048, 'Width',0.0062,'GroundPlaneLength', 0.0100, 'GroundPlaneWidth', 0.0100, 'PatchCenterOffset', [0 0],'FeedOffset', [0.0010 0],'Load', lumpedElement,'Tilt',[180],'TiltAxis',[0 1 0] );

H = zeros(NR,N);
cte = (4*pi)/lambda; 
cont=1;m=1;n=1;

for m = 1:NR
    for n=1:N
        vec_dif = AntenasRx(m,:)-RISv(n,:);
        [Azimut(m,n),Elevacio(m,n),Distancies(m,n)] = cart2sph(vec_dif(1),vec_dif(2),vec_dif(3));
        [eT, hTx] = EHfields(Antenna,f,vec_dif.');
        [eR, hRx] = EHfields(antRx,frx,(-vec_dif.'));
        fase2 = exp(1i*kl.*Distancies(m,n));%en contrafase per treure una de les dos de la polaritzacio
        H(m,n) = cte.*(Distancies(m,n)).* fase2.*(eT.'*eR);
    end
end
filename = string(D);
dAstr = string(dAnombre);
NRstr = string(NR);
filename = strcat('ARX_Cruzados_txd2M64x64_rxdRx12dmd_dRx16_D',filename);
save(filename,'H');
