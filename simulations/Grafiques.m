clear all
VecDis = 100*1.259.^[0:1:13];
VecDis1 = 100*1.259.^[0:1:13];
for Dist = 1:length(VecDis)
    valor_nom_corba = VecDis(Dist);
    nom_corba = string(valor_nom_corba);
    nom_corba = replace(nom_corba,'.',',');
    loadname = strcat("FAR_ARXO_Cruzados_A16_txd2M64x64_rxdRx12dmd_dRx16_P20,30_D",nom_corba);
    LH = load(loadname);
    MLP1(Dist) = LH.BER_MSE(2);
end

figure(1)
loglog(VecDis1,MLP1);grid; hold on;
xlabel('Distancia');
grid on;
ylabel('BER');
title('ML & MMSE 16 bits: RIS 64x64, RX arrays de 8 dipolos equiespaciados del centro de la RIS')
%legend('ML, far field, Amplitud a 1 m = 40','MSE, far field, Amplitud a 1 m = 40','ML, near field, Amplitud a 1 m = 40','MSE, near field, Amplitud a 1 m = 40')
legend('ML, SNR a 1m = 87dB','ML, SNR a 1m = 97dB',  ...
    'MMSE, SNR a 1m = 87dB','MMSE, SNR a 1m = 30dB')
%legend('Near field, Amplitud = 20', 'Near field, Amplitud = 30','Near field, Amplitud = 40','Far field, Amplitud = 20','Far field, Amplitud = 30','Far field, Amplitud = 40')
%legend('Near field, Amplitud = -10', 'Near field, Amplitud = 0','Near field, Amplitud = 10','Near field, Amplitud = 20','Near field, Amplitud = 30','Near field, Amplitud = 40')
