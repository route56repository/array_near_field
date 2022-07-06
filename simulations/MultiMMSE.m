
clear all;
%System parameters
Ntx = 64*64; %Transmit elements
Nrx = 4; %Receive antennas
Ls = 32*32;  %Elements per group
groups = Ntx/Ls; %Number of group elements

M = zeros(Ntx, groups); cont = 0; cont2 = 0;
for j = 1:sqrt(groups)   %definicio M
    for i = 1: sqrt(Ls)
        for m = 1: sqrt(groups)
            for n = 1:sqrt(Ls)
                M(n+cont*sqrt(Ls),m+cont2*sqrt(groups))=1;
            end
        cont = cont+1;
        end
    end
    cont2 = cont2+1;
end
%simulation parameters
N = 1e5;  %número de bits = groups*N

P = [-60];
VecDis = 100*1.259.^[0:1:13];
Autovalores = zeros(Nrx,length(VecDis));

for Dist = 1:length(VecDis) 
    valor_nom_corba = VecDis(Dist);
    nom_corba = string(valor_nom_corba);
    nom_corba = replace(nom_corba,'.',',');
    loadname = strcat("ARXO_SEPAMulti_A4_txd2M64x64dA_rxd10mdmd_D",nom_corba);
    %loadname = "";
    LH=load(loadname);
    H =LH.H;
    HM = H*M;
    A = HM'/(HM*HM');
    %BER computation
    for idx=1:length(P) % for each P (dB) value we compute the Pe
        a = 10^(P(idx)/20); %amplitud 
        s = sign(rand(groups,N)-0.5); %Symbol generation
        w = (randn(Nrx,N)+1i*randn(Nrx,N))/sqrt(2); %Noise
              
        m = sqrt(a^2/trace(A*A'));
        y_multi = HM*A*s*m+w;
        detectada_multi = sign(real(y_multi));
        err_multi = ne(detectada_multi,s);
        traza(Dist) = trace(A*A');
        BER(Dist,idx,1) =sum(sum(err_multi))/(groups*N);
        BER(Dist,idx,2) =sum(err_multi(1,:))/N;
        BER(Dist,idx,3) =sum(err_multi(2,:))/N;
        BER(Dist,idx,4) =sum(err_multi(3,:))/N;
        BER(Dist,idx,5) =sum(err_multi(4,:))/N;
    end
    Autovalores(:,Dist) = svd(H);
end
figure(1)
loglog(VecDis,BER(:,1,1));
grid; hold on;

xlabel('Distancia');
ylabel('BER');
title('ZF RIS 64x64 Rx 16 usuarios (dipolos) d_usuarios = 3m')
figure(2)
loglog(VecDis,traza);