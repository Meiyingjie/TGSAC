% Demo Program for the paper: 
%   SDRSAC: Semidefinite-Based Randomized Approach for Robust Point Cloud
%   Registration without Correspondences - Huu Le, Thanh-Toan Do, Tuan
%   Hoang and Ngai-Man Cheung - CVPR 2019
%       
%  
clc;clear;close all;
add_dependencies; 
% addpath(genpath('./RANSAC'));
dataset = 'synthetic';
method = 'SDRSAC';          % Without Correspondences
RR=[0.1125   -0.9588    0.2609
    0.8369    0.2330    0.4953
   -0.5357    0.1626    0.8286];
TT=[-0.12392735;-0.936076188;0.52748802];

% Read configuration containing hyperparameters
config = readConfig(dataset); 

% Read Data
load(config.matPath);
a=pcread(config.plyPath);
b=pcread(config.plyPathB);
MM=a.Location';
BB=b.Location';
flag=100000;
if length(MM)>=flag
    sap1=0.2;
else
    sap1=0;
end
if length(BB)>=flag
    sap2=0.2;
else
    sap2=0;
end
% sap1=0.15;
% sap2=0.15;
% figure;
% pcshow(a);
% hold on;
% pcshow(b);
% aa=pcdownsample(a,'gridAverage',0.8); 
% bb=pcdownsample(b,'gridAverage',0.8);  
aa=pcdownsample(a,'random',1); 
bb=pcdownsample(b,'random',0.4); 
% M=aa.Location';
M=load('C:\Users\23560\Desktop\point_model.txt')';
B=bb.Location';
% B=bb.Location'*1000;
% if length(MM)>=0.4e5
% aa1=pcdownsample(a,'gridAverage',0.12); 
% MM=aa1.Location';
% end
% if length(BB)>=0.4e5
% bb1=pcdownsample(b,'gridAverage',0.12); 
% BB=bb1.Location';

% end
% % 
% % ind=find(D(2,:)>0.08);
% % ind2=find(M(3,:)>0.03);
% % D=D(:,ind);
% % M=M(:,ind2);
% % D=D(:,1:100);
% % M=M(:,18000:20000);

% M=pcread('./data/ye.pcd');
% % M=double(M.Location)';
% % plot3(M(1,:), M(2,:), M(3,:), 'r.');
% % axis off;
% M=pcdownsample(M,'random',0.1);
% M=double(M.Location)';
% for i=1:size(M,1)
%     M(i,:)=(M(i,:)-repmat(min(M(i,:)), 1, size(M,2)))/(max(M(i,:))-min(M(i,:)))+1;
% %     M(i,:)=(M(i,:)-repmat(min(min(M)), 1, size(M,2)))/(max(max(M))-min(min(M)));
% %     M(i,:)=(M(i,:)-repmat(min(D(i,:)), 1, size(M,2)))/(max(D(i,:))-min(D(i,:)));
% end
%% 1ҶƬ
% Mm=RR*M + repmat(TT, 1, size(M,2));
% Dm=M';
% Mmm=Mm';
% noise=Mm+(randn(size(M,1),size(M,2))-0.5)*100/900;
% idx=randsample(size(noise,2),10000);
% B=[Mm,noise(:,idx)];
% MM=B';
% re_M=M';
% re_B=B';
%% 2ҶƬ
% Mm=RR*M + repmat(TT, 1, size(M,2));
% idx_Mm=randsample(size(Mm,2),500);
% Mm1=Mm(:,idx_Mm);
% Dm=M';
% Mmm=Mm';
% noise=Mm+(randn(size(Mm,1),size(Mm,2))-0.5)*100/900;
% idx=randsample(size(noise,2),9500);
% B=[Mm1,noise(:,idx)];
% MM=B';
% re_M=M';
% re_B=B';
%%
% ind=find(D(2,:)>1.6);
% ind2=find(Mm(2,:)>1.6);
% D=D(:,ind);
% M=Mm(:,ind2);

figure;
plotPointClouds(B, M, 'b.','r.')
axis equal;
% axis([-2,2,0,3,0,2])
% figure;
% plotPointClouds(B, Mm, 'b.','r.')

% pcshow(M');
% hold on;
% pcshow(B');
% figure;
% plot(M(1,:),M(2,:));
% Run SDRSAC
%%
% opt.method='affine'; % use rigid registration
% opt.viz=1;          % show every iteration
% opt.outliers=0.6;   % use 0.6 noise weight to add robustness 
% 
% opt.normalize=1;    % normalize to unit variance and zero mean before registering (default)
% opt.scale=1;        % estimate global scaling too (default)
% opt.rot=1;          % estimate strictly rotational matrix (default)
% opt.corresp=0;      % do not compute the correspondence vector at the end of registration (default)
% 
% opt.max_it=100;     % max number of iterations
% opt.tol=1e-8;       % tolerance
% 
% X=D';
% Y=M';
% X=X(5000:10000,:);
% Y=Y(5000:10000,:);
% % registering Y to X
% [Transform, Correspondence]=cpd_register(X,Y,opt);
% 
% figure,cpd_plot_iter(X, Y); title('Before');
% 
% % X(Correspondence,:) corresponds to Y
% figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');
%%
% bestRT=[ 9.88035971e-001    1.54222655e-001    5.41374166e-004    5.34167810e+000
% -1.54223348e-001    9.88034647e-001    1.64157616e-003   -2.56153095e+000
% -2.81728199e-004   -1.70542883e-003    9.99998506e-001   -1.46640894e+000];
% tran=bestRT(1:3,1:3)*M + bestRT(:,4);
out = pointCloudReg(B, M, B, M, config, method);
disp(out);
