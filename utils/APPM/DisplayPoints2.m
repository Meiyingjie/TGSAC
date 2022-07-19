function DisplayPoints2(Model, Scene);

set(gca,'FontSize',12);
plot3(Model(:,1),Model(:,2),Model(:,3),'k+');
hold on;
plot3(Scene(:,1),Scene(:,2),Scene(:,3),'go');
axis equal;
% xlim([min(Scene(:,1)),max(Scene(:,1))]);
% ylim([min(Scene(:,2)),max(Scene(:,2))]);   
% zlim([min(Scene(:,3)),max(Scene(:,3))]);   
pbaspect([1,1,1]);