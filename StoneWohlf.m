clear all;

dev=.00001;
dp=10^-3;

deg=45; theta=deg*pi/180;

 inu=1;
 for phiu=[dev*pi:dp*pi:pi-dev*pi]
  hh= -sin(2*(phiu-theta))./sin(phiu)/2;
     if cos(2*(phiu-theta))+hh*cos(phiu)>0
      angu(inu)=phiu;
      hu(inu)= -sin(2*(phiu-theta))./sin(phiu)/2;
      inu=inu+1;
     end
 end
 ind=1;
 for phid=pi+dev*pi:dp*pi:2*pi-dev*pi;
  hh= -sin(2*(phid-theta))./sin(phid)/2;
     if cos(2*(phid-theta))+hh*cos(phid)>0
      angd(ind)=phid;
      hd(ind)= -sin(2*(phid-theta))./sin(phid)/2;
      ind=ind+1;
     end
 end
 

mhu=cos(angu); mhu(1,length(mhu))=-1;
mhd=cos(angd); mhd(1,length(mhd))=1;

figure(1); 
subplot(2,1,1);plot(angu*180/pi,hu,angd*180/pi,hd);  
xlabel('\phi (deg)');ylabel('h');
axis([0 360 -2 2])
subplot(2,1,2);
plot(hu,mhu,hd,mhd); xlabel('h');ylabel('cos{\phi}');
axis([-1.1 1.1 -1.1 1.1])
