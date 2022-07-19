function [SM, SB, idx2] = getCorrespondences2(X,M,B)
 NN = size(M, 2); 
 pp=1:1:NN;
 rr=nchoosek(pp,3);
 %% 行表示M的索引，列表示B的索引
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
%% 利用tabulate函数统计点出现的频率，频率越高越可能是匹配点，选取频率最高的4个点
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