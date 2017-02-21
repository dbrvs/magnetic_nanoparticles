
r=115e-9; %distance between particle centers
mD=1e-24; %mass of spring
l=3e-9; %length spring
k=10*1.38e-23*300/pi^2; %spring const
m=2e-18; %mass of particle
R=56e-9; %particle radius

hbar=6.6e-34;

n=0:4;

Ev=sqrt(k/m).*(n+1/2);

Er=hbar*n.*(n+1)/(1/24*mD*l^2+1/4*m*r^2);

Et=hbar*n.*(n+1)/(4/5*m*R^2);

%plot(n,Ev,n,Er,n,Et)