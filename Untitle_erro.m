clc;clear;close all;
addpath(genpath('./utils'));
% a=pcread('data/bun045.ply');
% a=pcread('C:/Users/23560/Desktop/涡轮基mesh/mesh_ply.ply');
a=pcread('data/ye2.pcd');
bb=a.Location;
figure;
plot3(bb(:,1),bb(:,2),bb(:,3), '.','Color',[138, 43, 226]/255,'MarkerSize',20);
axis off
% aa=pcdownsample(a,'random',64/length(a.Location)); 
aa=pcdownsample(a,'random',400/length(a.Location)); 
MM=aa.Location;
tic;
tt=MyCrust(double(MM));
toc;
tri=delaunay(double(MM(:,1:3)));
shp=alphaShape(double(MM(:,1)),double(MM(:,2)),double(MM(:,3)));
% 
figure;
plot3(MM(:,1),MM(:,2),MM(:,3), '.','Color',[138, 43, 226]/255,'MarkerSize',20);
axis off
% plotPointCloud('data/bun045.ply');
figure;
trimesh(tt,double(MM(:,1)),double(MM(:,2)),double(MM(:,3)),'edgecolor',[138, 43, 226]/255,'LineWidth',1.5)
colormap autumn;
axis off
figure;
trisurf(tri,double(MM(:,1)),double(MM(:,2)),double(MM(:,3)))
colormap autumn;
% % figure;
% % trisurf(tri(:,1:3),double(MM(:,1)),double(MM(:,2)),double(MM(:,3)))
% % colormap autumn;

% clc;clear;close all
% gai=0.99;
% lv=0.01;

% lv=3/64;
% gai=1-lv;
% uyuy=0;
% m=3;
% n=3;
% % nchoosek(pp,3)
% for i=1:n
%     yuyu=vpa(gai^(m-i+1)*lv^(i-1)*nchoosek(m,i-1),10);
%     uyuy=yuyu+uyuy;
% %     yuyu=vpa(gai^(n-i+1)*lv^(i-1),10);
% %     uyuy=yuyu+uyuy;
% end
% uyuy
% log(0.01)/log(uyuy)
% lv*(gai^3+gai^2*lv+lv^2*gai)+gai^4


% factorial(16)/factorial(16-3)/factorial(3)
% base=importdata('./data/office/groundtruth/office.txt');
% tran=importdata('C:/Users/23560/Desktop/pointset/面积office点云坐标.txt');
% for i=1:8
%       R_base=base(1+(i-1)*4:3+(i-1)*4,1:3);
%       R_tran=tran(1+(i-1)*4:3+(i-1)*4,1:3);
%       T_base=base(1+(i-1)*4:3+(i-1)*4,4);
%       T_tran=tran(1+(i-1)*4:3+(i-1)*4,4);
% %       flag=(trace(R_base*R_tran')-1)/2;
% %       erro_R(i)=solve_R(flag);     
%       erro_R(i)=acosd((trace(R_base*R_tran')-1)/2);
%       erro_T(i)=norm(T_base-T_tran);
% end
% sum(erro_R)/8
% sum(erro_T)/8
% function pp=solve_R(flag)
% t=0.2*pi/180;
% pp=fsolve(@(x)myfun(x,flag),t);
% end
% % erro=abs(tran-base);
% % sum(erro)/64
% function F = myfun(x,flag)
% F = cos(x)-flag;
% end