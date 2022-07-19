% a=pcread('F:\新桌面\资料1\数字孪生\全新文献\点云\程序\数据\armadillo\ArmadilloBack_0.ply');
% b=pcread('F:\新桌面\资料1\数字孪生\全新文献\点云\程序\数据\armadillo\ArmadilloBack_60.ply');
% a=pcread('F:\新桌面\资料1\数字孪生\全新文献\点云\程序\数据\arch\s2.ply');
% % ttt=[2,3,4];
% % bbb=[3,6,7];
% % bbb-ttt
% % a=pcread('data\bun000.ply');
% b=pcread('F:\新桌面\资料1\数字孪生\全新文献\点云\程序\数据\arch\s4.ply');
% plyMaps('bunny10') = 'F:\2021年\数字孪生\盘模型\pan3pcd.pcd';
%     plyMapsB('bunny10') = 'F:\2021年\数字孪生\盘模型\real.ply';
a=pcread( 'F:\2021年\数字孪生\盘模型\pan3pcd.pcd');
b=pcread( 'F:\2021年\数字孪生\盘模型\real.ply');
MM=a.Location';
BB=b.Location';
% dsq = sum(power(MM(:,1:1000) - BB(:,1:1000), 2),1);
% ER=2/(mean(dsq)+2)
addpath('F:\新桌面\博士相关算法\SDRSAC-master - 副本\utils');
aa=pcdownsample(a,'random',1); 
bb=pcdownsample(b,'random',0.1); 
% aa=pcdownsample(a,'gridAverage',0.7); 
% bb=pcdownsample(b,'gridAverage',0.7);  
bb11=pcdownsample(b,'random',64/length(BB)); 
MM=bb.Location';
BB=load('C:\Users\23560\Desktop\point_model.txt')';
B_Tree = KDTreeSearcher(BB');
bestRT=[0.792867600917816,0.0207952484488487,-0.609040260314941,-0.209480747580528;-0.0913797318935394,0.992181956768036,-0.0850847586989403,0.875801980495453;0.602508127689362,0.123112447559834,0.788569271564484,1.77361917495728];
tran=bestRT(1:3,1:3)*MM + bestRT(:,4);
record=tran';
% save C:\Users\23560\Desktop\all.txt -ascii record;
inls_icp = countCorrespondences(tran, B_Tree, 0.015)
% params.algorithm = 'kdtree';
% params.trees = 8;
% params.checks = 64;
% % params.nareaPairThresh = 3;
% [srcIdx,dist] = flann_search(double(tran),double(BB),1,params);
% di=find(dist<=0.001);
% length(di)
% srcIdx
inls_icp/length(MM)
figure;
plotPointClouds(BB,tran, 'b.','r.')
% plot3(MM(1,:), MM(2,:), MM(3,:), 'r.', 'MarkerSize',0.8);
% axis off; 

% figure;
% mesh(D)
% max('*')
% en=class(':');
% xx='char'
% class(en)
% figure;
% pcshow(M')
% uu=[15,8,7,8,15,8,9];
% xx=tabulate(uu);
% [~,order]=sort(xx(:,3),'descend');
% roro=xx(:,1);
% roro2=roro(order);
% pp=1:1:15;
% rr=nchoosek(pp,3)
% size(rr)
% u1=[2 3 1;4 5 1;9 0 1];
% u2=[2 3 1;9 0 1;4 5 1];
% u3=[9 0 1;2 3 1;4 5 1];
% det(u1)
% det(u2)
% det(u3)
% 
% double count_triangle_area(Point a,Point b,Point c){
% 	double area = -1;
% 	 
% 	double side[3];//存储三条边的长度;
%  
% 	side[0] = sqrt(pow(a.x - b.x,2)+pow(a.y - b.y,2) + pow(a.z - b.z,2)); 
% 	side[1] = sqrt(pow(a.x - c.x,2)+pow(a.y - c.y,2) + pow(a.z - c.z,2));
% 	side[2] = sqrt(pow(c.x - b.x,2)+pow(c.y - b.y,2) + pow(c.z - b.z,2)); 
%  
% 	//不能构成三角形;
% 	if(side[0]+side[1]<=side[2] || side[0]+side[2]<=side[1] || side[1]+side[2]<=side[0]) return area; 
%  
% 	//利用海伦公式。s=sqr(p*(p-a)(p-b)(p-c)); 
% 	double p = (side[0]+side[1]+side[2])/2; //半周长;
% 	area = sqrt(p*(p-side[0])*(p-side[1])*(p-side[2])); 
% 	
% 	return area;
% }
% u2=[1 5;6 3];
% trace(u1'*u2*u1)
% trace(u1'*u1*u1)
% close all
% dataset = 'synthetic';
% method = 'SDRSAC';          % Without Correspondences
% 
% 
% % Read configuration containing hyperparameters
% config = readConfig(dataset); 4 
% 
% % Read Data
% load(config.matPath);
% a=pcread(config.plyPath);
% b=pcread(config.plyPathB);
% figure;
% pcshow(a);
% hold on;
% pcshow(b);
% bb=pcdownsample(b,'random',0.1);  
% aa=pcdownsample(a,'random',0.1);  
% 
% currCart=D';
% currCart=currCart(:,1:2);
% refCart=M';
% refCart=refCart(:,1:2);
% grid_num=10;
% % 这里需要把边界放宽一点
% % 这里需要把边界放宽一点
% x_up=    max([currCart(:,1);refCart(:,1)]);
% x_low=    min([currCart(:,1);refCart(:,1)]);
% y_up=    max([currCart(:,2);refCart(:,2)]);
% y_low=    min([currCart(:,2);refCart(:,2)]);
% gridSize=[(x_up-x_low)/(grid_num-1),(y_up-y_low)/(grid_num-1)];
% grid_dist_ref=cell(grid_num,grid_num);
% grid_dist_cur=cell(grid_num,grid_num);
% for j=1:size(refCart,1)
%     Pos=(refCart(j,:)-[x_low,y_low]);
%     Pos=[Pos(1)/gridSize(1),Pos(2)/gridSize(2)];
%     row=Pos(1)+1;
%     col=Pos(2)+1;
%     grid_dist_ref{ceil(row),ceil(col)}=[grid_dist_ref{ceil(row),ceil(col)};refCart(j,:)];
%     grid_dist_ref{ceil(row),floor(col)}=[grid_dist_ref{ceil(row),floor(col)};refCart(j,:)];
%     grid_dist_ref{floor(row),ceil(col)}=[grid_dist_ref{floor(row),ceil(col)};refCart(j,:)];
%     grid_dist_ref{floor(row),floor(col)}=[grid_dist_ref{floor(row),floor(col)};refCart(j,:)];
% end
% for j=1:size(currCart,1)
%     Pos=(currCart(j,:)-[x_low,y_low]);
%     Pos=[Pos(1)/gridSize(1),Pos(2)/gridSize(2)];
%     row=Pos(1)+1;
%     col=Pos(2)+1;
%     grid_dist_cur{ceil(row),ceil(col)}=[grid_dist_cur{ceil(row),ceil(col)};currCart(j,:)];
%     grid_dist_cur{ceil(row),floor(col)}=[grid_dist_cur{ceil(row),floor(col)};currCart(j,:)];
%     grid_dist_cur{floor(row),ceil(col)}=[grid_dist_cur{floor(row),ceil(col)};currCart(j,:)];
%     grid_dist_cur{floor(row),floor(col)}=[grid_dist_cur{floor(row),floor(col)};currCart(j,:)];
% end
% figure
% hold on
% for i=1:grid_num
%     for j=1:grid_num
%         Ref=grid_dist_ref{i,j};
%         Cur=grid_dist_cur{i,j};
%         if length( Ref)>0
%             plot(Ref(:,1),Ref(:,2),'*')
%         end
%         if length(Cur>0)
%             plot(Cur(:,1),Cur(:,2),'o')
%         end
%     end
% end
% figure,
% cellplot(grid_dist_ref)