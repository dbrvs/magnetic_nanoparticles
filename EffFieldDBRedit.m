function R53=EffFieldDBRedit(conc)
% FelderhofJones effective field approximation for magnetization 
% of magnetic nanoparticles from JBW/AW and edit by DBR2012

%ff=f;
tic;
%Initialize variables
Ho   = 20;   % field amp [mT]
f    = 10^3; % field freq [Hz]
visc = .001; % viscosity [Pa*s]
T    = 293;  % temperature [K] 
r    = 25;   % NP radius [nm] 5 or 57
Ms   = 70;   % NP saturation mag [emu/g] 31 or 70
tPts = 2^8;  % num time points
Thr  = .001; % ode threshold convergence
hmax = 10;   % number of harmonics
%Keff = 4.68e4;    %cubic anisotropy constant [J/m^3] GOYA


% constants
kB   = 1.38e-23; %boltzmann const. [J/K]
rho  = 3200;     %NP density [kg/m^3] from TDS 3.2g/ccm

%calculate variables
V     = 4/3*pi*(r*10^-9)^3; %NP volume [m^3]
mu    = Ms*V*rho;           %magnetic moment [J/T]
tau   = 3*visc*V/kB/T;      %Einstein relax time [s]
%tau   = 1e-10*exp(Keff*V/kB/T);

%conc=15; %[mg/mL]
ps=6e11*conc; %[particles/mL]
md=(10^-6/ps)^(1/3); %[m/particle]
c=mu/md^3;    % interaction depends on mu and mean distance


tspan = (0:tPts)/tPts/f;    %one point more than a cycle 
Ea    = Ho*mu/1000/kB/T;    %unitless field

%% Solve ODE

M=zeros(tPts+1,2); fM=zeros(tPts,2);
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

% get one first period
IC=1e-10; %initial condition very close to 0
[~,M(:,1)] =ode15s(@(t,m) SimMag(t,m,f,tau,Ea,c),tspan,IC,options);
fM(:,1)=fft(M(1:tPts,1));

% convergence loops checks between periods until meet tolerance
tol1=1; tol2=1; iter=2; MAXtol1=Thr; MAXtol2=Thr;
while tol1>MAXtol1 || tol2>MAXtol2;
   IC=M(tPts+1,1); %new initial condition is final point previous 
   
   % get next period
   [t,M(:,2)]=ode15s(@(t,m) SimMag(t,m,f,tau,Ea,c),tspan,IC,options);
   
   % compare periods and iterate
   fM(:,2)=fft(M(1:tPts,2));
   tol1=norm(sum(M(:,2)-M(:,1))./sum(M(:,1)));
   tol2=norm((fM(4,2)-fM(4,1))/fM(4,1));
   M(:,1)=M(:,2); fM(:,1)=fM(:,2);
   if tol1>MAXtol1 || tol2>MAXtol2; iter=iter+1; end;
   %dM=ifft(fM(:,2).*[0:size(fM,1)-1]');
end

%% Data analysis and plotting
pM=abs(fM(:,2)); H=zeros(1,hmax);
 for i =1:hmax; j=2*(i-1)+1; H(i)=pM(j+1)*j; end %find the odd derivative harmonic
R53=H(3)/H(2); %MSB signal

BB=sin(2*pi*f*t);
AA=abs(trapz(BB,M));
aa=AA(2);

%sar=mu*Ho*aa*f/rho*V;

%   %hold on;figure(3);figuresize(10,4,'inches' )
%    subplot(1,12,1:3); plot(t*1000,M)
%      xlabel('Time (ms)'); ylabel('Normalized mean magnetization')
%    subplot(1,12,5:6); plot(H)
%      xlabel('Harmonic #');  ylabel('Signal Magnitude');
%      xlim([0 hmax])
%    subplot(1,12,8:9); plot(BB,M)
%      xlabel('Applied Field (mT)'); ylabel('Normalized mean magnetization')
%      title(num2str(aa))
%    subplot(1,12,11:12); plot(f/10^6,aa,'x','MarkerSize',10); 
%      xlabel('Frequency (MHz)'); ylabel('Hysteresis loop area')
%           title(num2str(1/tau/10^6))
% 
% toc;

%m=M(:,2);


function dm=SimMag(t,m,f,tau,Ea,c) 
il = 3*m*(35-12*m^2)/(35-33*m^2); %Pade inverse Langevin
dm = -m/tau*(1-(Ea*sin(f*2*pi*t)+3*c*m)/il); %eff field approx

