clear all
tic
%System parameters
Ntx = 64*64; %Transmit elements
Nrx = 16; %Receive antennas
Ls = 16*16;  %Elements per group
groups = Ntx/Ls; %Number of group elements

%simulation parameters
N = 1e4;  %Number of symbols group to compute Pe (Joan las llama repeticiones)-> número de bits = groups*N
P = 100:-1:-100; % SNR at the receiver to test (in dBs)
Nruns = 1; % to experiment with different Rayleigh channels (channel is random: different channels may give different curves)

Pe = zeros(Nruns,length(P));
for experiment=1:Nruns %for each experiment we work with a different channel
    
    %Channels
    LH=load("txd2M64x64_rxdRx6dmd_dRx6_D10000"); %Rayleigh channel generation (unitary expected amplitude for each individual channel)
    %H=(randn(Nrx,Ntx)+1i*randn(Nrx,Ntx))/sqrt(2);
    %H = H/sqrt(mean(mean(abs(H)))); %si normalizamos hay menos cambios de canal a canal, pero ...
    ... es más correcto no hacerlo y aumentar Nruns para hacer la comparativa)
    H = LH.H*10000/1000;
    
      Heq2 = zeros(Nrx,groups); %Equivalent channel
      for idx=1:groups
         Heq2(:,idx)=sum(H(:,(idx-1)*Ls+1:idx*Ls),2);
      end
     
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
    HM = H*M;

    %Pe computation
    for idx=1:length(P) % for each P (dB) value we compute the Pe
        a = 10^(P(idx)/20); %amplitud 
        s = sign(rand(groups,N)-0.5); %Symbol generation
        x = kron(s,ones(Ls,1)); %Transmitted signal
        w = (randn(Nrx,N)+1i*randn(Nrx,N))/sqrt(2); %Noise
        y3 = H*a*x+w; %Received signal;
        y = HM*a*s+w; %para comprobar: y3 e y deberían ser iguales
        y2 = Heq2*a*s+w;
        %MMSE Receiver
        W = inv(a^2*(HM*HM')+eye(Nrx))*(HM*a);
        d = sign(real(W'*y)); %Detection
        err = ne(d,s);
        Pe(experiment,idx) =sum(sum(err))/groups/N;
    end
    
end
if Nruns>1
    Pe_all = mean(Pe);
else 
    Pe_all = Pe;
end
toc
semilogy(P,Pe_all);grid; hold on;
xlabel('SNR (dB) for single element-to-single antenna channel');
ylabel('BER');
title('MMSE. Bits 16 Media luna 64x64. Rx 16 dipolos')
svd = svd(HM)