%%% DBR 2013 -- LLG Model simulations %%%
function [tau,ff,bb,SAR,MM]=Neel(f,Bval)

ff=f;  
bb=Bval;

cycs = 3; %number of cycles until reach equil
  maxh   = 15; %number of harmonics to examine
%% Choose the physical conditions/variables %%
  pflg   = 1;         %plot flag
  %f      = 10^5;      %frequency [Hz]
  Bs     = 0;
  tPts   = 2^8;       %number time points
  N      = 1000;      %number of particles
  r      = 5e-9;   %avg radius [m]
  Ms     = 31;        %saturation magnetization of 5nm [emu/g] GOYA paper
  T      = 293;
  Keff   = 4.68e4;    %cubic anisotropy constant [J/m^3] GOYA
  al     = 1;
  gam    = 1.76086;   %gyromag ratio from Weiz [Hz/T]
  rr     = r*10;         %dipole field constant 
  ani    = [0 0 0];
  rho    = 3200;     %NP density [kg/m^3] from TDS 3.2g/ccm

 %calculate other constants based on the conditions above
  mu0  = 4*pi*10^-7;
  kT   = 1.38e-23*T;          %Boltzmann energy   
  V    = 4/3*pi*r^3;          %NP core volume [m^3]
  Ms   = Ms*rho;               %saturation magnetization [J/Tm^3] see pg108 notebook
  mu   = Ms*V;            %magnetic moment [J/T]
  Dn   = kT*al/mu/gam;
  K    = 2*Keff/Ms/mu0;
  C    = gam/(1+al^2);
  di   = mu/4/pi/rr^3;
  mass = rho*V;            % nanoparticle mass [kg]
  tau = 1e-10*exp(Keff*V/kT);

  
 m  = zeros(N,3); m(:,1) = 1; %make matrix initial mags. all in x
 %m = initrand(N); % random initial conditions see function for details
 B = zeros(N,3); B(:,3) = Bval/mu0/1000;%make matrix B in z dir.
 n = repmat(ani,N,1);
  
 per=1/f; dt=per/tPts; t=0:dt:per*cycs; %calculate timestep given cycles and timepoints
 in=1;  M=zeros(length(t),3); Heff=zeros(length(m),3);

 %% stochastic diffeq loop
for i=t;
  M(in,:)=mean(m);
  meanf=repmat(M(in,:),N,1); 

  h=normrnd(0,1,N,3); %stochastic term for fluctuations
  H=B*cos(2*pi*f*i); %drive field over time
    H(:,1)=Bs/1000/mu0;
  H2=B*cos(2*pi*f*(i+dt)); %drive field over time
    H2(:,1)=Bs/1000/mu0;

%    for jk=1:length(m);     
%     if norm(n)~=0;Heff(jk,:) = H(jk,:) + K*n*dot(n,m(jk,:));end
%    end

Heff(:,1) = H(:,1) + K*n(:,1).*m(:,1) + di*meanf(:,1);
Heff(:,2) = H(:,2) + K*n(:,2).*m(:,2) + di*meanf(:,2);
Heff(:,3) = H(:,3) + K*n(:,3).*m(:,3) + di*meanf(:,3);


   dmbar = C*(dt*(cross(m,Heff)+al*cross(m,cross(m,Heff))) +...
            sqrt(2*Dn*dt)*(cross(m,h)+al*cross(m,cross(m,h))));
   mbar = m - dmbar;
   dm = C/2*(dt*(cross(m,Heff)+al*cross(m,cross(m,Heff))...
         +cross(mbar,H2)+al*cross(mbar,cross(mbar,H2)))+...
          sqrt(2*Dn*dt)*(cross(m,h)+al*cross(m,cross(m,h))+...
           cross(mbar,h)+al*cross(mbar,cross(mbar,h))));
   m=m-dm; 
   %normalize the magnitude  
 nm = sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2); 
 m(:,1)=m(:,1)./nm; m(:,2)=m(:,2)./nm; m(:,3)=m(:,3)./nm;
 
 in=in+1;
end;
  
%% analysis/plot
  lastper=2*tPts+1:3*tPts; harm=zeros(maxh,3);
   BB=cos(2*pi*lastper*dt*f); MM=M(lastper,3);
   AA=abs(trapz(BB,MM)); SAR=AA*mu*Bval/mass/1000;
  pM=abs(fft(M(lastper,:)));
   for i =1:10; harm(i,:)=pM(i+1,:)*i; end %find the derivative harmonics
  %R53=harm(5,3)/harm(3,3); %MSB signal
  %R42=harm(4,1)/harm(2,1); %sMSB signal
 
   if pflg ~= 0
  %magnetization and harmonics plot
  figure(1); %figuresize(8,4,'inches' )
   subplot(2,2,1:2); plot(t*1000,M)
     xlabel('Time (ms)'); ylabel('Normalized mean M')
     legend('M_x','M_y','M_z')
   subplot(2,2,3); bar(harm)
     xlabel('Harmonic #');  ylabel('Signal Magnitude');
     xlim([0 maxh+1]); set(gca, 'XTick', 0:5:maxh);
   subplot(2,2,4)%hysteresis plot
   plot(BB,MM)
   xlabel('Normalized Applied Field'); 
   ylabel('Norm Mean M')
    title(num2str(SAR))
  end
