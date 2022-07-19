function X=solveSDP3(W,config)
X=zeros(size(W,1),size(W,2));
for i=1:config.nareaPairThresh
    [m1,idx]=max(W);
    % idy=1:size(a,2);
    [~,idx2]=max(m1);
    % a(idx(idx2),idx2)
    X(idx(idx2),idx2)=1;
    W(idx(idx2),:)=-1000;
    W(:,idx2)=-1000;
end
