%% sweep fields and freqs
function [mz,r,F,A]=sweepB

nb=2;
nf=5;

m0=[1 0 0];
chT=10^-4;
pflg = 1;
cycs = 3;      %number of cycles
N    = 10^4;   %number of particles
tPts = 2^8;   %time points

mz=zeros(tPts,nf,nb);
r=zeros(nf,nb);

for b=1:nb;   
    for f=1:nf;
        F(f)=f*500; A(b)=9+b; a0=[0 0 A(b)];

        [mz(:,f,b),r(f,b)]=THB(a0,F(f),chT,m0,cycs,N,tPts); 
    end
end

%%
%save('data',mz,r)
%S=ScaleCorr(F,r(:,1),r(:,2));
YETR=sqrt(1+.21*A(1).^2)./sqrt(1+.21*A(2).^2);
MYT1=2*atanh(max(mz(:,:,1)))/A(1);
MYT2=2*atanh(max(mz(:,:,2)))/A(2);
MYTR=MYT2./MYT1;

plot(F,r(:,1),F,r(:,2),F*YETR,r(:,2),F.*MYTR,r(:,2))
legend('\alpha_1','\alpha_2','scaled \alpha_2')
xlabel('Frequency (Hz)'); ylabel('Harmonic ratio')
%title(['Scale value = 'num2str(S)])

%%% DBR 2013 -- Theory Brownian code %%%
function [mz,R53]=THB(alpha,freq,tau,mi,cycs,N,tPts)

 m  = repmat(mi,N,1);    %matrix of initial mags
 a  = repmat(alpha,N,1); %matrix of fields
 dt = 1/freq/tPts;  
 tf = cycs/freq; 
 M  = zeros(tPts,3);

in=1;
for i=0:dt:tf
  %mm(in,:,:)=m; %for single particle trajectories
  M(in,:)=mean(m);  
  A=a*cos(2*pi*freq*i); %drive field over time
  
 m=m+(cross(cross(m,A),m)-2*m)*dt/tau+...
      cross(randn(N,3),m)*sqrt(dt/tau);
  
%normalize the magnitude  
 nm = sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2); 
 m(:,1)=m(:,1)./nm; m(:,2)=m(:,2)./nm; m(:,3)=m(:,3)./nm;
  
in=in+1;
end
  lastper=2*tPts+1:3*tPts; H=zeros(10,3);
  mm=M(lastper,:);
  mz=mm(:,3);
  pM=abs(fft(mm));
  for i =1:10; H(i,:)=pM(i+1,:)*i; end %find the derivative harmonics
  R53=H(5,3)/H(3,3); %MSB signal
  %R42=H(4,1)/H(2,1); %sMSB signal
  %hx=H(:,1);   
 
    