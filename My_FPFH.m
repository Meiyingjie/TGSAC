function FPFH=My_FPFH(P,nor_p,idx_fe_p,r,idx_rn_p,dis_p)
% P 原始点云数据 3×n
% nor_p 原始点云的法向量矩阵 3×n
% idx_fe_p 点云特征点在原数据中的位置索引向量
% r FPFH邻域参数
% [idx_rn_p,dis_p]为rangesearch函数的输出参数，可缺省
% dis_rneighbor 对应r邻域点到中心点的欧氏距离
if nargin < 4
    error('no bandwidth specified')
end
if nargin < 5
    Mdl_p = createns(P','NSMethod','kdtree','Distance','minkowski','p',2);
%     size(Mdl_p)
    [idx_rn_p,dis_p]=rangesearch(Mdl_p,P',r);
end

FPFH = zeros(length(idx_fe_p),33);
loop = 0;
for i=idx_fe_p
    loop = loop + 1;
    idx_rneighbor = idx_rn_p{i};
    dis_rneighbor = dis_p{i};
    SPFH = My_SPFH(P,nor_p,idx_rneighbor,dis_rneighbor);
    SPFH2 = zeros(1,33);
    for j = 2:length(idx_rneighbor)
        idx_rneighbor2 = idx_rn_p{idx_rneighbor(j)};
        dis_rneighbor2 = dis_p{idx_rneighbor(j)};
        rnei_SPFH = My_SPFH(P,nor_p,idx_rneighbor2,dis_rneighbor2);
        weight = dis_rneighbor(j);
        SPFH2 = SPFH2+rnei_SPFH./weight;
    end
    FPFH(loop,:) = SPFH + SPFH2./(length(idx_rneighbor)-1);
end
end


function SPFH=My_SPFH(P,nor_p,idx_rneighbor,dis_rneighbor)
% P 原始点云数据
% nor_p 原始点云的法向量矩阵
% idx_rneighbor 本次r邻域点在原始点云中的位置索引向量，其中第一个元素为邻域中心点的位置索引
% dis_rneighbor 对应r邻域点到中心点的欧氏距离

rneighbor = P(:,idx_rneighbor);
nor_rneighbor = nor_p(:,idx_rneighbor);
u = nor_rneighbor(:,1);
SPFH = zeros(1,33);
for j=2:length(idx_rneighbor)
    v = cross(u,(rneighbor(:,j)-rneighbor(:,1))./dis_rneighbor(j));
    w = cross(u,v);
    alpha = dot(v,nor_rneighbor(:,j));
    vote = uint8(ceil((alpha+1)*(11/2)));
    if vote<=0
        vote=1;
    elseif vote>11
        vote = 11;
    end
    SPFH(vote)=SPFH(vote)+1;
    phi = dot(u,(rneighbor(:,j)-rneighbor(:,1))./dis_rneighbor(j));
    vote = uint8(ceil((phi+1)*(11/2)));
    if vote<=0
        vote=1;
    elseif vote>11
        vote = 11;
    end
    SPFH(vote+11)=SPFH(vote+11)+1;
    theta = sin(atan(dot(w,nor_rneighbor(:,j))./dot(u,nor_rneighbor(:,j))));
    vote = uint8(ceil((theta+1)*(11/pi)));
    if vote<=0
        vote=1;
    elseif vote>11
        vote = 11;
    end
    vote = uint8(vote);
    SPFH(vote+22)=SPFH(vote+22)+1;
end
end