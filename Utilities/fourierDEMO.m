%%
x=0:.001:2*pi;

for i=1:15;
    
y(i,:)=exp(-i+1)*sin(i*x);

end

yy = y([1 2 5],:);
sumy= yy(1,:) + yy(2,:) + yy(3,:);

figure(1); figuresize(8,4,'inches' )
   subplot(1,7,1:3); plot(x,yy)
     xlabel('x'); ylabel('f(x)')
     xlim([0 2*pi])
     %legend('Mx','My','Mz','Location','NorthEast');
     legend('boxoff')
   subplot(1,7,5:7); plot(x,sumy)
     xlabel('x'); ylabel('\Sigma f(x)')
     xlim([0 2*pi])

     fancyGraph(gcf)