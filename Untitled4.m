pc= pcread('C:\Users\23560\Desktop\pointset\new\Project1\Project1\data\out.pcd');
pp=pc.Location;
pp=reshape(pp,size(pp,1)*size(pp,2),[]);

P=pp';
P=P(:,150000:200000);

[index1,index11,~]=find(isnan(P));
P(:,index11)=[];
nor_pp=pc.Normal;
nor_pp=reshape(nor_pp,size(nor_pp,1)*size(nor_pp,2),[]);
nor_p=nor_pp';
nor_p=nor_p(:,150000:200000);
nor_p(:,index11)=[];
[index2,index22,~]=find(isnan(nor_p));
nor_p(:,index22)=0;
figure;
pcshow(P');
% FPFH=My_FPFH(P,nor_p,1:64,2);