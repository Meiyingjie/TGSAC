clc
clear
% close all
% P = ascread('bun000.asc');
% P=P{2};
Pp =pcread('F:\新桌面\资料1\数字孪生\全新文献\点云\程序\数据\facade\s2.ply');
MM=Pp.Location';
aa=pcdownsample(Pp,'gridAverage',0.2); 
% MM2=aa.Location';
% aa=pcdownsample(aa,'random',10000/length(MM2)); 
P = aa.Location';
figure;
x=P(1,1:1:end);
y=P(2,(1:1:end));
z=P(3,(1:1:end));
c=z+1;
scatter3(x,y,z,1,c,'filled');
colorbar
view(2)
title('原始数据')
Mdl_p = createns(P','NSMethod','kdtree','Distance','minkowski','p',2);
[idx_rn_p,dis_p]=rangesearch(Mdl_p,P',0.5);
idx_fe_p = My_ISS(P,0.5 ,0.8,0.4,idx_rn_p,dis_p);


figure;
x=P(1,idx_fe_p(1:1:end));
y=P(2,idx_fe_p(1:1:end));
z=P(3,idx_fe_p(1:1:end));
c=z+1;

scatter3(x,y,z,2,c,'filled');
colorbar
view(2)
title('关键点选取结果')
