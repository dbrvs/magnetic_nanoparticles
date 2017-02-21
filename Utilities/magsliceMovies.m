function magsliceMovies(mmx,mmz)


imax=10;
tPts=length(mmz);

figure(2); figuresize(8.5,5,'inches')
subplot(2,imax,1:imax)%,'Position',[.1 .8 .8 .1])
t=(0:tPts-1)/tPts;
B=cos(2*pi*t);

Mx=mean(mmx,2); Mz=mean(mmz,2);
plot(t,Mx,t,Mz,t,B)
xlim([0 max(t)])
xlabel('Fraction of cycle')
ylabel('Mean normalized M')
legend('M_x','M_x','B_~')
legend('Location','SouthEast')
legend('boxoff')


subplot(2,imax,1+imax)%,'Position',[xx yy ww hh]);
subplot(2,imax,5+imax)%,'Position',[xx yy ww hh]);
ylabel('M_z');

for i=1:imax
    n=i*floor(tPts/imax);


    if i==2;ylabel('M_x'); %set(gca, 'YTick', [-1 0 1]); 
    end
    if i==imax/2+1; xlabel('M_z');end;

subplot(2,imax,i+imax)%,'Position',[i/10 .1 .3 .3]);

scatter(mmx(n,:),mmz(n,:),4,'k');

%title([num2str(n/f*1000) ' ms']);

set(gca, 'XTick', []);
set(gca, 'YTick', []); 
box on
xlim([-1 1])
ylim([-1 1])

hold on; plot([-1 1], [0 0], 'k-'); 
plot([0 0], [-1 1], 'k-'); hold off;
end

fancyGraph(gcf)

%%
% figure(2)
% figuresize(3,3 ,'inches' )
% for j=1:length(mmz)
%     
% scatter(mmx(j,:),mmz(j,:),4,'blue');
% set(gca, 'XTick', [-1 0 1]);
% set(gca, 'YTick', [-1 0 1]); 
% box on
% xlim([-1 1])
% ylim([-1 1])
% xlabel('M_x/|M|')
% ylabel('M_z/|M|')
% 
% hold on; plot([-1 1], [0 0], 'k-'); 
% plot([0 0], [-1 1], 'k-'); hold off;
% 
% F(j) = getframe(gcf);
% end

%%

