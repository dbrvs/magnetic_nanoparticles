% Binding simulation 2013 DBR

function [mz,R53] = BIND(concsp)

tic;
f    = 400;   %frequency [Hz]
Bv   = 10;     %oscillating field amp [mT]
%K    = 10^-19;%binding order of pN*nm [Nm]
tPts = 2^8;    %time points
T    = 293;    %temperature [K]
rhy  = 57e-9;  %hydrodynamic radius [m]
rco  = 25e-9;  %core radius [m]
rho  = 3200;   %density [kg/m^3] from data sheet 3.2g/ccm
N    = 10^4;   %number of particles
visc = .001025;   %viscosity [Pas]
Ms   = 68;    %saturation magnetization [emu/g]
Bs   = 0;      %static field [mT]
%Bc   = 0;      %circular field component [mT]

  kT  = 1.38e-23*T;   %Boltzmann energy   
  Vhy = 4/3*pi*rhy^3; %NP hydrodynamic volume [m^3]
  Vco = 4/3*pi*rco^3; %NP core volume [m^3]
  mu  = rho*Ms*Vco;   %magnetic moment [J/T]
  gm  = 6*visc*Vhy;   %drag coefficient [Nsm]
  tE  = gm/2/kT;      %Einstein relaxation time [s]
   
  K=kT*200;
  
  Dv=kT/gm; Du=mu/gm; Dk=K/gm; %variables

 m  = initrand(N); %init mags randomly
 %m   = zeros(N,3); m(:,1) = 1; %init mags in x
 %meq = initrand(N); %random springs
 %meq = zeros(N,3); meq(:,2) = 1; %directional springs
 
 numsp=round(N*concsp);
 meq = zeros(N,3); meq(1:numsp,:) = initrand(numsp); %some random springs
 
 b   = zeros(N,3); b(:,3)=Bv/1000; %make matrix b in z dir.
 
 %calculate timesteps
 if f==0; dt=10^-5; tf=5*tE; 
 else per=1/f; dt=per/tPts; tf=3*per; 
 end
 t=0:dt:tf; 
 
%% Stochastic differential equation loop
in=1;    M=zeros(length(t),3); TB=zeros(N,3); TS=TB; Ts=TB;
for i= 0:dt:tf
  %mmx(in,:)=m(:,1); mmy(in,:)=m(:,2); mmz(in,:)=m(:,3); %for slices

    M(in,:)=mean(m);
  h=randn(N,3); %stochastic term for fluctuations
  B=b*cos(2*pi*f*i); %drive field over time
    B(:,1)=Bs/1000;
    %B(:,2)=Bc/1000*sin(2*pi*f*i); %circular field

  TB=Du*cross(cross(m,B),m); %mag torque
  
  if K~=0
  kk=Dk*cross(cross(m,meq),m); %spring torque direction
   
  %spring torque magnitude linear with angle
  angbtw=acos(dot(m,meq,2)); %angle from equilibrium
  
  TS(:,1)=kk(:,1).*angbtw; 
  TS(:,2)=kk(:,2).*angbtw; 
  TS(:,3)=kk(:,3).*angbtw;
  
  %for a=1:N; TS(a,:)=kk(a,:)*acos(dot(m(a,:),meq(a,:),2)./norm(m(a,:))); end
  
  TB=TB+TS;
  end
  
  Ts=sqrt(2*Dv)*cross(h,m); %stochastic torque
    
%stochastic diffeq Stratanovich with drift term
  dm = (TB-2*m*Dv)*dt + Ts*sqrt(dt);
  m = m + dm; 
%normalize
 nm = sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2); 
 m(:,1)=m(:,1)./nm; m(:,2)=m(:,2)./nm; m(:,3)=m(:,3)./nm;
 
 in=in+1;
end;

%% Plotting and harmonic analysis
  lastper=2*tPts+1:3*tPts; H=zeros(10,3);
  pMx=abs(fft(M(lastper,1)));
  pMy=abs(fft(M(lastper,2)));
  pMz=abs(fft(M(lastper,3)));%/length(M);
   for i =1:10; 
       H(i,1)=pMx(i+1)*i;
       H(i,2)=pMy(i+1)*i;
       H(i,3)=pMz(i+1)*i;
    end %find the derivative harmonic
    mz=M(:,3);
  R53=H(5,3)/H(3,3); %MSB signal
  %R42=H(4,1)/H(2,1); %sMSB signal
  
  figure;figuresize(8,4,'inches' )
   subplot(1,7,1:4); plot(t*1000,M)
     xlabel('Time (ms)'); ylabel('Normalized mean magnetization')
     legend('Mx','My','Mz','Location','NorthWest');
   subplot(1,7,6:7); bar(H)
     xlabel('Harmonic #');  ylabel('Signal Magnitude');
     legend('Hx','Hy','Hz','Location','NorthWest'); 
     %legend('boxoff')
toc;




  