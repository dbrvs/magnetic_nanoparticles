function vectormov(M)


figure()
%[x, y, z] = sphere(30);
%mesh(x,y,z);colormap(gray(120))
%set(gca, 'XTick', []);
%set(gca, 'YTick', []); 
%box on;
hold on;
grid on;
view(135,20)
% xlim([min(min(M)) max(max(M))])
% ylim([min(min(M)) max(max(M))])
% zlim([min(min(M)) max(max(M))])

xlim([-1 1]); ylim([-1 1]); zlim([-1 1])

xlabel('M_x'); ylabel('M_y'); zlabel('M_z')


axis('square')


%%
for i=1:length(M)
    

    plot3(M(i,1),M(i,2),M(i,3),'.','MarkerSize',8)

%F(i) = getframe(gcf);

end


