function [R,X,B] = twoloops(ni,a,ni2,a2)


%find the magnetic field of loops of wire

%help from I. Teliban, modified by D. Reeves 6/21/11
clf; %clears the current figure

%vaccuum permeability [Vs/Am, H/m, N/A^2, Tm/A, Wb/(A·m)]
mu0=.39*4*pi*10^-5; %-5 for cm .39 for inches

%viewing window

xmin=-1.5;xmax=1.5;
rmin=-2;rmax=2;
step=.01;

%set up the sizes for matrices in advance. doesn't actually seem to save
%that much time though... comment out for now
%xscale=(xmax-xmin+1)/step;rscale=(rmax-rmin+1)/step;

%B=zeros(xscale,rscale);Br=zeros(xscale,rscale);Bx=zeros(xscale,rscale);
%Br2=zeros(xscale,rscale);Bx2=zeros(xscale,rscale);X=zeros(xscale,rscale);
%R=zeros(xscale,rscale);X2=zeros(xscale,rscale);R2=zeros(xscale,rscale);

%first loop variables
%ni=2000;
%a=2;
xshift=-.75;
rshift=0;

%second loop variables
%ni2=0;
%a2=2;
x2shift=.75;
r2shift=0;

%calculate first loop
B0=(ni.*mu0)./(2.*a);
in1=1;
for x=(xmin):step:(xmax)
    in2=1;
    for r=(rmin):step:(rmax)
        if (r~=0)
            al=abs((r)./a);
            be=abs((x+xshift)./a);
            ga=((x+xshift)./r);
            q=((1+al).^2+be.^2);
            k=sqrt(4*al./q);
            [K,E] = ellipke(k.^2);
            Bx(in1,in2)=B0.*(1./(pi.*sqrt(q))).*(E.*(1-al.^2-be.^2)./(q-4.*al)+K);
            Br(in1,in2)=B0.*(ga./(pi.*sqrt(q))).*(E.*(1+al.^2+be.^2)./(q-4.*al)-K);            
        else
            Bx(in1,in2) = B0.*(a./sqrt(a.^2+x.^2)).^3;
            Br(in1,in2) = 0;
   
        end
        X(in1,in2)=x;
        R(in1,in2)=r;
        in2=in2+1;
    end
    in1=in1+1;
end

%second loop calculations
%want to shift the ring back by a certain amount so add this to x,
%this makes it seem farther away, field is weaker

B02=(ni2.*mu0)./(2*a2);
in3=1;
for x2=xmin:step:xmax%keep same array size as first loop so Matlab is happy
    in4=1;
    for r2=rmin:step:rmax
        if (r2~=0)
            al2=abs((r2+r2shift)./a2);
            be2=abs((x2+x2shift)./a2);
            ga2=((x2+x2shift)./r2);
            q2=((1+al2).^2+be2.^2);
            k2=sqrt(4*al2./q2);
            [K2,E2] = ellipke(k2.^2);
            Bx2(in3,in4)=B02.*(1./(pi*sqrt(q2))).*(E2.*(1-al2.^2-be2.^2)./(q2-4.*al2)+K2);
            Br2(in3,in4)=B02.*(ga2./(pi*sqrt(q2))).*(E2.*(1+al2.^2+be2.^2)./(q2-4.*al2)-K2);            
        else
            Bx2(in3,in4)=B02.*(a2./sqrt(a2.^2+x2.^2)).^3;
            Br2(in3,in4)=0;
        end
        X2(in3,in4)=x2;%keep same values so matrices match up except values are shifted in both B's respectively
        R2(in3,in4)=r2;
        in4=in4+1;
    end
    in3=in3+1;
end

%turn coil to face in the other direction
%Bx2=Bx2';
%Br2=Br2';
%Bx=Bx';
%Br=Br';

Bx=Bx2+Bx;
Br=Br2+Br;

B=sqrt(Bx.^2+Br.^2);
%surfc(X,R,B);
%% 
contour(R,X,B,200);
shading flat;
colormap('jet');

hidden off;
colorbar;
hold on;
%quiver(X,R,Bx,Br,5,'k');
%h=streamslice(X,R,Bx,Br);
%set(h,'Color','w')
xlabel('x (in)');
ylabel('y (in)');
title('Magnetic Field (T)');
caxis auto;

end
