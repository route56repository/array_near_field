clear all;
%System parameters
Ntx = 64*64; %Transmit elements
Nrx = 16; %Receive antennas
Ls = 16*16;  %Elements per group
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
N = 1e6;  %número de bits = groups*N

P = [105,107,110];

VecDis = 10*1.259.^[0:1:13];
Autovalores = zeros(Nrx,length(VecDis));
    
LH=load("MultiCurvado_txd2M64x64_Rxdsind(3)_dmd16_Dgrupos100_400.mat");
H =LH.H;
HM = H*M;
A = HM/(HM.'*HM);
BER = zeros(groups,length(P));
%BER computation
for idx=1:length(P) % for each P (dB) value we compute the Pe
    a = 10^(P(idx)/20); %amplitud 
    s = sign(rand(groups,N)-0.5); %Symbol generation
    w = (randn(Nrx,N)+1i*randn(Nrx,N))/sqrt(2); %Noise
                
    m = sqrt(a^2/trace(A*A'));
    y_multi = HM.'*A*m*s+w;
    detectada_multi = sign(real(y_multi));
    err_multi = ne(detectada_multi,s);
    
    BER(:,idx) =sum(err_multi.')/(N);
    BER(Dist,idx,2) =sum(err_multi(1,:))/N;
    BER(Dist,idx,3) =sum(err_multi(2,:))/N;
    BER(Dist,idx,4) =sum(err_multi(3,:))/N;
    BER(Dist,idx,5) =sum(err_multi(4,:))/N;

    traza = trace(A*A');
end

figure(1)
usuarios = linspace(1,16,16);
semilogy(usuarios,BER(:,1));grid; hold on;
semilogy(usuarios,BER(:,2))
semilogy(usuarios,BER(:,3))
xlabel('Distancia');
ylabel('BER');
title('ZF RIS 64x64 Rx 16 usuarios (dipolos) d_usuarios = 10m')

figure(2)
loglog(VecDis,traza);


