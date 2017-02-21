function M= detNeel
tic;

al=.01;
cycs=4;
f=1000;

m=[1; 0; 0];
n=[1; 0; 0]; %if norm(n)~=0; n=n/norm(n); end
B=[0; 0; f*10];

per=cycs/f;
dt=per/10000;

tspan=0:dt:per-dt;

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
%%

[~,M] =ode15s(@(t,m) SimMag(t,m,al,B,n,f),tspan,m,options);

%  in=1;
% for i=0:dt:cycs
%  M(in,:)=m;
%   H=B*cos(2*pi*i)+n*dot(n,m);
%   m=m-(cross(m,H)+al*cross(cross(m,H),m))*dt/(1+al^2);
%  in=in+1;
% end

%lastper=(cycs-1)/dt:cycs/dt;
%%
t=0:dt:per-dt;
BB=cos(2*pi*f*t);
%SAR=trapz(M(lastper,3));

figure(1); 
subplot(2,1,1);
plot(t,M);
xlabel('Time');ylabel('Magnetization');
subplot(2,1,2);
plot(BB,M,'r');
xlabel('Field');ylabel('Magnetization');
%title(['SAR = ' num2str(SAR)]);

toc

function dm=SimMag(t,m,al,B,n,f) 
   H=B*cos(2*pi*f*t)+n*dot(n,m);

dm=-cross(m,H)-al*cross(m,cross(m,H));
