function idx_feature = My_ISS(p, r, e1,e2,idx,dis)
if nargin < 4
    error('no bandwidth specified')
end
if nargin < 5
    Mdl = createns(p','NSMethod','kdtree','Distance','minkowski','p',2);
    [idx,dis] = rangesearch(Mdl,p',r);
end
numpts = size(p,2);
flag = zeros(1,numpts);
% tzz = zeros(size(p));
for i = 1:numpts
    if length(idx{i})<2
        continue
    end
    x = p(:,idx{i}(2:end));%r邻域点坐标
    w = 1./dis{i}(2:end);
    p_bar = p(:,i);% 中心点坐标
    P = repmat(w,3,1).*(x - repmat(p_bar,1,size(x,2))) * ...
        transpose(x - repmat(p_bar,1,size(x,2))); %spd matrix P
    P = P./sum(w);
%     if any(isnan(P(:)))
%         save debug.mat 
%     end
    [~,D] = eig(P);
    lam = sort(abs(diag(D)),'descend');% 三个特征值由大到小排列
    if lam(2)/lam(1)<=e1 &&lam(3)/lam(2)<e2
        flag(i)=1;
    end
%     tzz(:,i)=lam;
end
% tzz(1,:)=tzz(2,:)./tzz(1,:);
% tzz(2,:)=tzz(3,:)./tzz(2,:);
% tzz(3,:)=[];
idx_feature = find(flag);
end