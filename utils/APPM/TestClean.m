clc;clear;
% M = rand(30,2)*5; 
% S = TransformPoint([1,1,.25].*randn(1,3),M); 
% a=pcread('cloud_bin_0.ply');
% b=pcread('cloud_bin_1.ply');
% bb=pcdownsample(b,'random',0.1);  
% aa=pcdownsample(a,'random',0.1);  
% D1=bb.Location;
% D2=aa.Location;
load('bunny10.mat')
D1=M';
D2=D';
% D1=D1(:,1:2);
% D2=D2(:,1:2);
figure(3); clf;
KCReg2(D1,D2,2,1);
