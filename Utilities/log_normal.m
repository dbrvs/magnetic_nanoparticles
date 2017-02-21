
mu=1;
sigma=.01;
N=1000000;

x=exp(mu+sigma^2/2)

R = lognrnd(mu,sigma,N,1);

%hist(R,200)

r=60e-9*x;%mean(R);

V=4/3*pi*r^3;

mu=76*25*V;