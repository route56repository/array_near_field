%Diferentes distancias
NR = 16;
N = 64*64;
agrupacio = 16*16;
Ls = 16;
agrupacio_lateral = sqrt(agrupacio);
sectors = N/agrupacio;
sectors_lateral = sqrt(sectors);

%en distancia fija
P = [-40,-20];
repeticiones = 1e4; %1e6
%Mat_Hs = (randn(NR,N)+1i*randn(NR,N))/sqrt(2);
%Mat_Hs = load('');

Mat_R = zeros(NR,2^sectors);
seq_recibida_ML = zeros(sectors, repeticiones);
Errors_ML = zeros(sectors,repeticiones,200);
seq_recibida_MSE = zeros(sectors, repeticiones);
Errors_MSE = zeros(sectors, repeticiones,200);
BER_MSE = zeros(length(P),1);
BER_ML = zeros(length(P),1);
M = zeros(N, sectors); cont = 0; cont2 = 0;
for j = 1:sectors_lateral   %definicio M
    for i = 1: agrupacio_lateral
        for m = 1: sectors_lateral
            for n = 1:agrupacio_lateral
                M(n+cont*agrupacio_lateral,m+cont2*sectors_lateral)=1;
            end
        cont = cont+1;
        end
    end
    cont2 = cont2+1;
end

sequencia_candidatas = zeros(sectors,1);
Mat_seq_cand = zeros(sectors, 2^sectors);
cont =1;
for n = 1:2^sectors %generacio candidates
    Mat_seq_cand(:,n) = sequencia_candidatas(:);
    if cont < 2^sectors
        for m = 1:log2(cont)+1
            sequencia_candidatas(m) = bitget(cont,m);
        end
    end
    cont = cont+1;
end


cont = 0;
Distmax = 1000;
Resdis = 100;
ResPot = 1;
SNR = 40; %dB
experiments = 1;

VecDis = 100*1.259.^[0:1:13];



for Dist = 1:length(VecDis)
valor_nom_corba = VecDis(Dist);
nom_corba = string(valor_nom_corba);
nom_corba = replace(nom_corba,'.',',');
loadname = strcat("ARXO_txd2M64x64dA_rxd12dmd_D",nom_corba,".mat");
Mat_Hs = load(loadname);
H = Mat_Hs.H;

    for idx =1: length(P)
        %Variacio Potencia
        amp = 10^(P(idx)/20);
        HM = H*M*amp;
        Mat_seq_rx = HM*(2*Mat_seq_cand-1);
        cont = 1;
        for rep =1:repeticiones
            
            seq_enviada = randi([0 1], sectors,1);
            sequencia = seq_enviada*2-1;
            w = (randn(NR,1)+1i*randn(NR,1))/sqrt(2);
            r = HM*sequencia+w;
            
            %ML
            Mat_R = repmat(r,1,2^sectors);
            ML = vecnorm(Mat_R - Mat_seq_rx).^2;
            [MLmod, MLrx] = min(ML);
            seq_recibida_ML(:,cont) = (de2bi(MLrx-1,sectors)); %,false int2bit
            Errors_ML(:,cont,idx) = ne(seq_recibida_ML(:,cont),seq_enviada);

            %MSE
            
            A = ((HM*HM')+eye(NR))\(HM);
            MSE = A'*r;
            seq_recibida_MSE(:,cont) = (sign(real(MSE))+1)/2;
            Errors_MSE(:,cont,idx) = ne(seq_recibida_MSE(:,cont),seq_enviada);

            cont = cont+1;
        end

    Total_Errors_ML = sum(Errors_ML(:,:,idx), 'all');
    Total_Errors_MSE = sum(Errors_MSE(:,:,idx), 'all');
    BER_ML(idx) = Total_Errors_ML/(sectors*repeticiones);
    BER_MSE(idx) = Total_Errors_MSE/(sectors*repeticiones);
    %valuePot = string(P(idx));
    %valuePot = erase(valuePot,".");
    %sexp = string(exp);
    %resname = strcat('sepA16Q_txMd2dA100_rxM16_ML_MSE_A16_','P_',valuePot, '_exp',sexp);
    %save(resname,"BER_ML","BER_MSE");
    %resname = strcat('FAR_ArrayLineal_txMd2_rxdmd16_ML_MSE_A16_D',nom_corba);
    %save(resname,"BER_ML","BER_MSE");
    end
resname = strcat('ARXO_QSEPA16_txd2M64x64dA100_rxd12dmd_P-40,-20_D',nom_corba);
save(resname,"BER_ML","BER_MSE");    
end