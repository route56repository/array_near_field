

VecDis = 100*1.259.^[0:1:13];
LHfar = load("FARCruzados_txd2M64x64_rxdRx12dmd_dRx6_D100000");
for Dist = 1:length(VecDis)
    valor_nom_corba = VecDis(Dist);
    nom_corba = string(valor_nom_corba);
    nom_corba = replace(nom_corba,'.',',');
    loadname = strcat("Cruzados_txd2M64x64_rxdRx12dmd_dRx6_D",nom_corba);
    LHnear = load(loadname);
    Hfar = LHfar.H*1e5/VecDis(Dist);
    Garray(Dist) = norm(LHnear.H)/norm(Hfar);
end

semilogx(VecDis,Garray);grid; hold on;
xlabel('Distancia');
grid on;
ylabel('Nromalaized Array Gain');
title("Array Gain Near field / Far Field")