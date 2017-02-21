clear all;

dev=.00001;
dp=10^-3;

 deg1=45; t1=deg1*pi/180;
 deg1=45; t2=deg1*pi/180;

 
Ms=1;
B=1;
V=1;
K=1;
d=1;

eb=Ms*V*B;
ek=K*V;
ed=10^-7*Ms^2*V^2/d^3;

h=zeros(1/dp,1/dp);
p2=[dev*pi:dp*pi:pi-dev*pi];
in1=1;
for p1=[dev*pi:dp*pi:pi-dev*pi];
 

%for equilibrium 

bb=sin(p1)-sin(p2);
kk=ek/eb*(sin(2*(p1-t1))-sin(2*(p2-t2)));
dd=ed/eb*(3*cos(p1).*sin(p2)+sin(p1-p2)-...
    3*cos(p2).*sin(p1)+sin(p1-p2));

h(in1,:)=(kk+dd)./bb;

in1=in1+1;    
end







