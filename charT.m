
%%
V=4/3*pi*(56e-9)^3;
Vcore=4/3*pi*(25e-9)^3;
visc=.001;
%B=[.000001;2;5;8;10;15;25]./1000;
B=.0001:.001:.05;
mu=3200*75*Vcore;
kB=1.38e-23;
T=293;

eT=3*visc*V/kB/T*10^6;

%alpha0 = mu*B/kB/T;
alpha0=0:.1:20;

msat=coth(alpha0)-1./alpha0;

chT=2*atanh(msat)./alpha0;

%chT=log((1+msat)./(1-msat))./alpha0;

YET=1./sqrt(1+.21*alpha0.^2);

plot(alpha0,chT,'-',alpha0,YET,'--')
xlabel('\alpha')
ylabel('\tau/\tau_B')
