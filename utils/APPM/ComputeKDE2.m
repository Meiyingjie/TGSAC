%compute the Kernel Density Estimation defined by the N by 2 array P
function KDE = ComputeKDE2(P)

global resolution;
global min_val;
global max_val;

%find the range of grid points
grids = round( (max_val-min_val)/resolution)+20;
grids1=grids(:,1:2);
KDE1 = zeros(grids1);
KDE2 = zeros(grids1);
KDE3 = zeros(grids1);
% size(grids1)
% size(grids1)
% size(KDE1)
% kernel = zeros(grids(:,:,1));
% size(grids)
start = min_val - 10*resolution*ones(1,3);

%find the points whose centers are within range
index = find( P(:,1)>=min_val(1)- 6*resolution & P(:,1) <= max_val(1)+ 6*resolution & ...
    P(:,2)>=min_val(2)- 6*resolution & P(:,2)<=max_val(2)+ 6*resolution & ...
    P(:,3)>=min_val(3)- 6*resolution & P(:,3)<=max_val(3)+ 6*resolution);

for ii = 1:length(index)
    point = P(index(ii),:);
    center = round( (point-start)/resolution);
    x_range = center(1)-3:center(1)+3;
    y_range = center(2)-3:center(2)+3;
    z_range = center(3)-3:center(3)+3;
%     size(KDE1)
%     size(KDE2)
%     size(KDE3)
    if max(x_range)>=min(size(KDE1,1),size(KDE1,2))
        [index1,index11,~]=find(x_range>=min(size(KDE1,1),size(KDE1,2)));
        x_range(index1,index11)=min(size(KDE1,1),size(KDE1,2));
    end
    if max(y_range)>=min(size(KDE1,1),size(KDE1,2))
        [index2,index22,~]=find(y_range>=min(size(KDE1,1),size(KDE1,2)));
        y_range(index2,index22)=min(size(KDE1,1),size(KDE1,2));
    end
    if max(z_range)>=min(size(KDE1,1),size(KDE1,2))
        [index3,index33,~]=find(z_range>=min(size(KDE1,1),size(KDE1,2)));
        z_range(index3,index33)=min(size(KDE1,1),size(KDE1,2));
    end
%     max(x_range)
%     max(y_range)
%     max(z_range)
    x_val = start(1) + x_range * resolution - point(1);
    y_val = start(2) + y_range * resolution - point(2);
    z_val = start(3) + z_range * resolution - point(3);
    kernel_x = exp(- x_val.*x_val/(resolution*resolution));
    kernel_x = kernel_x/sum(kernel_x);
    kernel_y = exp(- y_val.*y_val/(resolution*resolution));
    kernel_y = kernel_y/sum(kernel_y);
    kernel_z = exp(- z_val.*z_val/(resolution*resolution));
    kernel_z = kernel_z/sum(kernel_z);
%     size(kernel_x)
%     size(KDE)
%     for flag = 1:max(size(kernel_z,1),size(kernel_z,2))
%         kernel(1:max(size(kernel_z,1),size(kernel_z,2)),1:max(size(kernel_z,1),size(kernel_z,2)),flag)=(kernel_x'*kernel_y)*kernel_z(:,flag);
%     end
% %     size(kernel)
%     KDE(x_range,y_range,z_range) = KDE(x_range,y_range,z_range) + kernel(flag,flag,flag);
%     KDE1
%     KDE2
%     KDE3
%     x_range
%     y_range
    KDE1(x_range,y_range) = KDE1(x_range,y_range) + kernel_x'*kernel_y;
    KDE2(x_range,z_range) = KDE2(x_range,z_range) + kernel_x'*kernel_z;
    KDE3(y_range,z_range) = KDE3(y_range,z_range) + kernel_y'*kernel_z;
end

nm1 = sqrt(sum(sum(KDE1.^2)));
nm2 = sqrt(sum(sum(KDE2.^2)));
nm3 = sqrt(sum(sum(KDE3.^2)));
KDE1 = KDE1/nm1;
KDE2 = KDE2/nm2;
KDE3 = KDE3/nm3;
KDE=[KDE1 KDE2 KDE3];