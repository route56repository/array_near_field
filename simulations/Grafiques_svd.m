clear all
VecDis = 100*1.259.^[0:1:13];
Nrx = 4;
%LHfar = load("Cruzados_txd2M64x64_rxdRx12dmd_dRx6_D10000");
%FARd1e5=load("FARCruzados_txd2M64x64_rxdRx12dmd_dRx6_D100000");
lhfar = load("FAR1e5_ARXO_Cruzados_txd2M64x64_rxdRx12dmd_dRx16_D100000");
normfar = norm(lhfar.H*1e5);

for Dist = 1:length(VecDis)
    valor_nom_corba = VecDis(Dist);
    nom_corba = string(valor_nom_corba);
    nom_corba = replace(nom_corba,'.',',');
     loadname = strcat("ARXO_Cruzados_txd2M64x64_rxdRx12dmd_dRx16_D",nom_corba);
     LH = load(loadname);
     near_svd(Dist,:) = svd(LH.H);
     far_svd(Dist,:) = svd(lhfar.H*100000/VecDis(Dist));
     %far_svd2(Dist,:) = svd(FARd1e5.H*100000/VecDis(Dist));
%      Hfar = lhfar.H*1e5/VecDis(Dist);
%      far_svd(Dist,:) = svd(Hfar);
     normnear(Dist) = norm(LH.H);
      normfarv(Dist) = normfar/VecDis(Dist);
end
c = 299792458;
lambda = c/30e9;

%PropH = FARd1e5.H*(1e5/1e4)./LHfar.H;
%absPropH = abs(PropH);
%fasesH = angle(PropH);
relacio_svdnear = 20*log10(near_svd(:,10)./near_svd(:,1));
relacio_svdfar = 20*log10(far_svd(:,6)./far_svd(:,1));
%relacio_svdfar2 = 20*log10(far_svd2(:,5)./far_svd2(:,1));
faseCalculada = angle(exp(-1i*2*pi*90000/lambda));

%relacion = relacio_svdfar2./relacio_svdfar;

plotvec = VecDis.';
x = linspace(1,Nrx,Nrx);
% 
% figure(1);
% %semilogx(VecDis, relacio_svdnear);
% semilogx(VecDis, relacio_svdfar);
% %semilogx(VecDis, relacio_svdfar2);
% grid; hold on;
% xlabel('Distancia');
% grid on;
% ylabel('relación autovalores (dB)');
% title('relación autovalores: tx RIS, rx 2 arrays de 8 dipolos')
% legend('lambda2/lambda1','lambda4/lambda1','lambda6/lambda1','lambda4/lambda1','lambda5/lambda1')

%figure(2)
%semilogx(VecDis, near_svd);

figure(3)
semilogx(VecDis,normnear./normfarv)
title('Relación normas en campos near field/far field. Receptor array linial de 16 dipolos')
xlabel('Distancia');
ylabel('Relación normas near field/far field')
grid on

figure(4)
norm_svd = (far_svd(:,1)).^2/(far_svd(1,1))^2;
semilogx(VecDis,norm_svd);
xlabel('Distancia');
ylabel('Normalaized array gain')
title('Array Gain')
grid on 
hold on
