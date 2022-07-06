%Diferentes distancias
NR = 16;
N = 64*64;
agrupacio = 16*16;
Ls = 16;
agrupacio_lateral = sqrt(agrupacio);
sectors = N/agrupacio;
sectors_lateral = sqrt(sectors);

P = [20];
repeticiones = 1e3; 

Mat_R = zeros(NR,2^sectors);
seq_recibida_MSE = zeros(sectors, repeticiones);
Errors_MSE = zeros(sectors, repeticiones,200);
BER_MSE = zeros(length(P),1);

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
        for rep =1:repeticiones
            seq_enviada = randi([0 1], sectors,1);
            sequencia = seq_enviada*2-1;
            w = (randn(NR,1)+1i*randn(NR,1))/sqrt(2);
            r = HM*sequencia+w;
            %MSE
            A = ((HM*HM')+eye(NR))\(HM);
            MSE = A'*r;
            seq_recibida_MSE(:,rep) = (sign(real(MSE))+1)/2;
            Errors_MSE(:,rep,idx) = ne(seq_recibida_MSE(:,rep),seq_enviada);   
        end
    Total_Errors_MSE = sum(Errors_MSE(:,:,idx), 'all');
    BER_MSE(idx) = Total_Errors_MSE/(sectors*repeticiones);
    %valuePot = string(P(idx));
    %valuePot = erase(valuePot,".");
    %sexp = string(exp);
    %resname = strcat('sepA16Q_txMd2dA100_rxM16_ML_MSE_A16_','P_',valuePot, '_exp',sexp);
    %save(resname,"BER_ML","BER_MSE");
    %resname = strcat('FAR_ArrayLineal_txMd2_rxdmd16_ML_MSE_A16_D',nom_corba);
    %save(resname,"BER_ML","BER_MSE");
    end
resname = strcat('ARXO_MSE_QSEPA16_txd2M64x64dA100_rxd12dmd_P20_D',nom_corba);
save(resname,"BER_MSE");    
end