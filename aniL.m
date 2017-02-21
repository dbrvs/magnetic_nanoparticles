function aniL
for i=1:100

    B(i)=.0002*i;%10^(i-1);

    aM(i) = calcAM(B(i));


end
plot(B,aM)


function aM = calcAM(B)

mu=4.5e-18;
K=4.68e-4; %J/m^3
kT=293*1.38e-23;
V=4/3*pi*(30e-9)^3; %core volume

a=mu*B/kT;
b=K*V/kT;

mi=(a-2*b)/2/sqrt(b);
pl=(a+2*b)/2/sqrt(b);

MI=mfun('dawson',mi);
PL=mfun('dawson',pl);


aM = (a*MI-a*exp(2*a)*PL+(exp(2*a)-1)*sqrt(b))/b/(exp(2*a)*PL-MI);

