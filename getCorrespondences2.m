function [SM, SB, idx2] = getCorrespondences2(X,M,B)
 NN = size(M, 2); 
 pp=1:1:NN;
 rr=nchoosek(pp,3);
 %% �б�ʾM���������б�ʾB������
 [w,d]=find(X==1);
 p=[];
 q=[];
 for i=1:length(w)
     p=[p rr(w(i),:)];
     q=[q rr(d(i),:)];
 end
%  for ii=1:length(d)
%      q=[q rr(d(ii),:)];
%  end
%% ����tabulate����ͳ�Ƶ���ֵ�Ƶ�ʣ�Ƶ��Խ��Խ������ƥ��㣬ѡȡƵ����ߵ�4����
% xx1=tabulate(p);
% xx2=tabulate(q);
% [~,order1]=sort(xx1(:,3),'descend');
% [~,order2]=sort(xx2(:,3),'descend');
% roro1=xx1(:,1);
% roro2=xx2(:,1);
% roro11=roro1(order1);
% roro22=roro2(order2);
% idx1=roro11(1:4);
% idx2=roro22(1:4);
[idx1,idx2]=TAB(p,q,3);
SM=M(:,idx1);
SB=B(:,idx2);   
end