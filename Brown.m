%%% DBR 2013 -- Torque balance simulations %%%
function [mz]=Brown(Bv,f)

pflg = 0;
cycs = 5;      %number of cycles
N    = 10^4;   %number of particles
tPts = 2^10;   %time points

T    = 293;    %temperature [K]
%f    = 10;   %frequency [Hz]
%Bv   = 10;   %oscillating field amp [mT]
Bs   = 0;       %static field [mT]
visc = .001;    %viscosity [Pas]
rhy  = 57e-9;   %hydrodynamic radius [m]
rco  = 25e-9;   %core radius [m]
rho  = 3200;    %density [kg/m^3] from data sheet 3.2g/ccm
Ms   = 70;      %saturation magnetization [emu/g]=10^-3[J/T]

  kT  = 1.38e-23*T;   %Boltzmann energy   
  Vhy = 4/3*pi*rhy^3; %NP hydrodynamic volume [m^3]
  Vco = 4/3*pi*rco^3; %NP core volume [m^3]
  mu  = rho*Ms*Vco;   %magnetic moment [J/T]
  gm  = 6*visc*Vhy;   %drag coefficient [Nsm]
  tE  = gm/2/kT;      %Einstein relaxation time [s]
    
  Dv=kT/gm; Du=mu/gm; %variables

%% variables for characteristic time
  alpha0 = mu*Bv/kT/1000;
  msat   = coth(alpha0)-1/alpha0;
  charT  = 2*atanh(msat)/alpha0*tE*10^6; %[µs]

%% initialize magnetizations, fields and times
 %m  = initrand(0,1,N); %init mags randomly
 m   = zeros(N,3); m(:,1) = 1; %init mags in x
 b   = zeros(N,3); b(:,3)=Bv/1000; %make matrix b in z dir.
 
 %calculate timesteps
 if f==0;  tf=tE/2; dt=tf/tPts;
 else per=1/f; dt=per/tPts; tf=cycs*per; 
 end
 t=0:dt:tf; t=10^6*reshape(t,length(t),1);
 
%% Stochastic differential equation loop
in=1;    M=zeros(length(t),3); %TB=zeros(N,3); Ts=TB;
for i=0:dt:tf
    
    %mm(in,:,:)=m; %for single particle trajectories
  
  M(in,:)=mean(m);  
  h=randn(N,3); %stochastic term for fluctuations
  B=b*cos(2*pi*f*i); %drive field over time
    B(:,1)=Bs/1000;
    %B(:,2)=Bc/1000*sin(2*pi*f*i); %circular field

  TB=Du*cross(cross(m,B),m); %mag torque  
  Ts=sqrt(2*Dv)*cross(h,m); %stochastic torque
    
%stochastic diffeq Stratanovich with drift term
  dm = (TB-2*m*Dv)*dt + Ts*sqrt(dt);
  m = m + dm; 

%normalize the magnitude  
 nm = sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2); 
 m(:,1)=m(:,1)./nm; m(:,2)=m(:,2)./nm; m(:,3)=m(:,3)./nm;
  
  in=in+1;
end;

%% Plotting and harmonic analysis
   %Mmax=max(M(:,3));

if f~=0;
  lastper=2*tPts+1:3*tPts; H=zeros(10,3);
  mm=M(lastper,:);
  mz=mm(:,3);
  pM=abs(fft(mm));
  for i =1:10; H(i,:)=pM(i+1,:)*i; end %find the derivative harmonics
  R53=H(5,3)/H(3,3); %MSB signal
  %R42=H(4,1)/H(2,1); %sMSB signal
  %hx=H(:,1);    
    
end

  
  if pflg ~= 0
  figure(1); figuresize(8,4,'inches' )
   subplot(1,7,1:4); plot(t*1000,M)
     xlabel('Time (ms)'); ylabel('Normalized mean magnetization')
     legend('Mx','My','Mz','Location','NorthEast');
     legend('boxoff')
   subplot(1,7,6:7); bar(H)
     xlabel('Harmonic #');  ylabel('Signal Magnitude');
     legend('Hx','Hy','Hz','Location','NorthEast'); 
     legend('boxoff')
     xlim([0 10])
  end



  