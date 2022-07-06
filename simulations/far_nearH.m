%10000vs100
c = 299792458;
lambda = c/30e9;

Hfar64 = load("D10000c.mat");
Hnear64 = load("D100000c.mat");
PropH64 = Hnear64.H*(100000/10000)./Hfar64.H;
absH64 = abs(PropH64);
fases = angle(PropH64);
faseCalculada = angle(exp(-1i*2*pi*90000/lambda));
save("ProporcionFases",'fases','faseCalculada');

Hfar5 = load("D10000_Rx5_N5.mat");
Hnear5 = load("D300_Rx5_N25.mat");
PropH5 = Hfar5.H.*(10000/300)./Hnear5.H;
absH5 = abs(PropH5);
absHnear5 = abs(Hnear5.H);
absHfar5 = abs(Hfar5.H);