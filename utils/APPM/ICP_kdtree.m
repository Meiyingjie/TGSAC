function [TR, TT, er, t] = ICP_kdtree(q,p,config,varargin)
% Perform the Iterative Closest Point algorithm on three dimensional point
% clouds.
%
% [TR, TT] = icp(q,p)   returns the rotation matrix TR and translation
% vector TT that minimizes the distances from (TR * p + TT) to q.
% p is a 3xm matrix and q is a 3xn matrix.
%
% [TR, TT] = icp(q,p,k)   forces the algorithm to make k iterations
% exactly. The default is 10 iterations.
%
% [TR, TT, ER] = icp(q,p,k)   also returns the RMS of errors for k
% iterations in a (k+1)x1 vector. ER(0) is the initial error.
%
% [TR, TT, ER, t] = icp(q,p,k)   also returns the calculation times per
% iteration in a (k+1)x1 vector. t(0) is the time consumed for preprocessing.
%
% Additional settings may be provided in a parameter list:
%
% Boundary
%       {[]} | 1x? vector
%       If EdgeRejection is set, a vector can be provided that indexes into
%       q and specifies which points of q are on the boundary.
%
% EdgeRejection
%       {false} | true
%       If EdgeRejection is true, point matches to edge vertices of q are
%       ignored. Requires that boundary points of q are specified using
%       Boundary or that a triangulation matrix for q is provided.
%
%
% Matching
%       {bruteForce} | Delaunay | kDtree
%       Specifies how point matching should be done. 
%       bruteForce is usually the slowest and kDtree is the fastest.
%       Note that the kDtree option is depends on the Statistics Toolbox
%       v. 7.3 or higher.
%
% Minimize
%       {point} | plane | lmaPoint
%       Defines whether point to point or point to plane minimization
%       should be performed. point is based on the SVD approach and is
%       usually the fastest. plane will often yield higher accuracy. It 
%       uses linearized angles and requires surface normals for all points 
%       in q. Calculation of surface normals requires substantial pre
%       proccessing.
%       The option lmaPoint does point to point minimization using the non
%       linear least squares Levenberg Marquardt algorithm. Results are
%       generally the same as in points, but computation time may differ.
%
% Normals
%       {[]} | n x 3 matrix
%       A matrix of normals for the n points in q might be provided.
%       Normals of q are used for point to plane minimization.
%       Else normals will be found through a PCA of the 4 nearest
%       neighbors.
%
% ReturnAll
%       {false} | true
%       Determines whether R and T should be returned for all iterations
%       or only for the last one. If this option is set to true, R will be
%       a 3x3x(k+1) matrix and T will be a 3x1x(k+1) matrix.
%
% Triangulation
%       {[]} | ? x 3 matrix
%       A triangulation matrix for the points in q can be provided,
%       enabling EdgeRejection. The elements should index into q, defining
%       point triples that act together as triangles.
%
%
% Weight
%       {@(match)ones(1,m)} | Function handle
%       For point or plane minimization, a function handle to a weighting 
%       function can be provided. The weighting function will be called 
%       with one argument, a 1xm vector that specifies point pairs by 
%       indexing into q. The weighting function should return a 1xm vector 
%       of weights for every point pair.
%
% WorstRejection
%       {0} | scalar in ]0; 1[
%       Reject a given percentage of the worst point pairs, based on their
%       Euclidean distance.
%
% Martin Kjer and Jakob Wilm, Technical University of Denmark, 2012

% Use the inputParser class to validate input arguments.
inp = inputParser;

inp.addRequired('q', @(x)isreal(x) && size(x,1) == 3);
inp.addRequired('p', @(x)isreal(x) && size(x,1) == 3);

inp.addOptional('iter', 10, @(x)x > 0 && x < 10^5);

inp.addParameter('Boundary', [], @(x)size(x,1) == 1);

inp.addParameter('EdgeRejection', false, @(x)islogical(x));

validMatching = {'bruteForce','Delaunay','kDtree'};
inp.addParameter('Matching', 'bruteForce', @(x)any(strcmpi(x,validMatching)));

validMinimize = {'point','plane','lmapoint'};
inp.addParameter('Minimize', 'point', @(x)any(strcmpi(x,validMinimize)));

inp.addParameter('Normals', [], @(x)isreal(x) && size(x,1) == 3);

inp.addParameter('NormalsData', [], @(x)isreal(x) && size(x,1) == 3);

inp.addParameter('ReturnAll', false, @(x)islogical(x));

inp.addParameter('Triangulation', [], @(x)isreal(x) && size(x,2) == 3);

inp.addParameter('Verbose', false, @(x)islogical(x));

inp.addParameter('Weight', @(x)ones(1,length(x)), @(x)isa(x,'function_handle'));

inp.addParameter('WorstRejection', 0, @(x)isscalar(x) && x > 0 && x < 1);

inp.addParameter('translation', false,  @(x)islogical(x));
inp.parse(q,p,varargin{:});
arg = inp.Results;
clear('inp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual implementation

% Allocate vector for RMS of errors in every iteration.
global Scene;
global Model;
global display_it;
t = zeros(arg.iter+1,1); 
%% �����趨
a=[2,1/2,1/4,0,-1/4,-1/2,-1,-2,-4,-8,-16,-32,-64,-128,-inf];%% �˻����
a_flag=1;
% Start timer
% tic;
% Transformed data point cloud
pt = p;
k=0;
% if length(p)<length(q)
%     p_flag=q;
%     q_flag=p;
% else
%     p_flag=p;
%     q_flag=q;
% end
Np = min(size(p,2),size(q,2));
% Allocate vector for RMS of errors in every iteration.
ER = zeros(arg.iter+1,1); 
stop=0;
% Initialize temporary transform vector and matrix.
T = zeros(3,1);
R = eye(3,3);
% Initialize total transform vector(s) and rotation matric(es).
TT = zeros(3,1, arg.iter+1);
TR = repmat(eye(3,3), [1,1, arg.iter+1]);
if length(p)<=length(q)
    if length(p)>length(q)
       p_flag=q;
       q_flag=pt;
    else
       p_flag=pt;
       q_flag=q;
    end
    % If Minimize == 'plane', normals are needed
    if (strcmp(arg.Minimize, 'plane') && isempty(arg.Normals))
        arg.Normals = lsqnormest(q_flag,4);
    end
    % If Matching == 'Delaunay', a triangulation is needed
    if strcmp(arg.Matching, 'Delaunay')
        DT = delaunayTriangulation(transpose(q_flag));
    end
    % If Matching == 'kDtree', a kD tree should be built (req. Stat. TB >= 7.3)
    if strcmp(arg.Matching, 'kDtree')
        kdOBJ = KDTreeSearcher(transpose(q_flag));
    end
end
%%
%disp('Running ICP ....');
% Go into main iteration loop
while k<=arg.iter && stop==0
    k=k+1;   
    if length(p)>length(q)
       p_flag=q;
       q_flag=pt;
    else
       p_flag=pt;
       q_flag=q;
    end
    %% ������ʽѡ��
    if length(p)>length(q)
        % If Minimize == 'plane', normals are needed
        if (strcmp(arg.Minimize, 'plane') && isempty(arg.Normals))
            arg.Normals = lsqnormest(q_flag,4);
        end
        % If Matching == 'Delaunay', a triangulation is needed
        if strcmp(arg.Matching, 'Delaunay')
            DT = delaunayTriangulation(transpose(q_flag));
        end
        % If Matching == 'kDtree', a kD tree should be built (req. Stat. TB >= 7.3)
        if strcmp(arg.Matching, 'kDtree')
            kdOBJ = KDTreeSearcher(transpose(q_flag));
        end
    end
    %% 
    % If edge vertices should be rejected, find edge vertices
    if arg.EdgeRejection
        if isempty(arg.Boundary)
            bdr = find_bound(q, arg.Triangulation);
        else
            bdr = arg.Boundary;
        end
    end

    % Do matching
    switch arg.Matching
        case 'bruteForce'
            [match,mindist] = match_bruteForce(q_flag,p_flag);
        case 'Delaunay'
            [match,mindist] = match_Delaunay(q_flag,p_flag,DT);
        case 'kDtree'
            [match,mindist] = match_kDtree(q_flag,p_flag,kdOBJ);
%             [match mindist] = match_kDtree2(q,pt,config);
    end

    % If matches to edge vertices should be rejected
    if arg.EdgeRejection
        p_idx = not(ismember(match, bdr));
        q_idx = match(p_idx);
        mindist = mindist(p_idx);
    else
%         p_idx = false(1, Np);
        p_idx = true(1, Np);
        q_idx = match;
    end
    
    % If worst matches should be rejected
    if arg.WorstRejection
        edge = round((1-arg.WorstRejection)*sum(p_idx));
        pairs = find(p_idx);
        [~, idx] = sort(mindist);
%         fifi=mindist;
%         max(mindist)
%         min(mindist)
        p_idx(pairs(idx(edge:end))) = false;
        q_idx = match(p_idx);
        mindist = mindist(p_idx);
%         fifi(p_idx)=10;
%         max(fifi)
%         min(fifi)
    end
%     size(q_idx)
%     edge
    if k == 1
        ER(k) = sqrt(sum(mindist.^2)/length(mindist));
    end
    if length(p)>length(q)
        q_idx_flag=p_flag(:,p_idx)';
        p_idx_flag=q_flag(:,q_idx)';
    else
        q_idx_flag=q_flag(:,q_idx)';
        p_idx_flag=p_flag(:,p_idx)';
    end
    Model = q_idx_flag;
    Scene = p_idx_flag;
    switch arg.Minimize
        case 'point'
            % Determine weight vector
            weights = arg.Weight(match);
%             weights1 =match;
             weights1 =mindist;
            [R,T] = eq_point(Model',Scene', weights(p_idx));
        case 'plane'
            weights = arg.Weight(match);
            [R,T] = eq_plane(Model',Scene',arg.Normals(:,q_idx),weights(p_idx));
        case 'lmapoint'
            weights1 =match;
            [R,T] = eq_lmaPoint(Model',Scene');
    end
    
    % Add to the total transformation
    TR(:,:,k+1) = R*TR(:,:,k);
    TT(:,:,k+1) = R*TT(:,:,k)+T;
    if arg.translation
        k_flag=80;
    else
        k_flag=15;
    end
    if k<=k_flag
        TR(:,:,k+1)=eye(3,3);
%         T=eq_point_T(q(:,q_idx),pt(:,p_idx), weights(p_idx));
    end
    % Apply last transformation
    pt = TR(:,:,k+1) * p + repmat(TT(:,:,k+1), 1, size(p,2));

    if mod((k-1),3)==0
        a_mark=a(a_flag);
        a_flag=a_flag+1;
        if a_flag>length(a)
            a_flag=length(a);
        end
    end
%     ER(k+1)=rms_error2(q(:,q_idx), pt(:,p_idx),a_mark,weights1(p_idx));
    ER(k+1)=rms_error2(Model',Scene',-inf,weights1);
    
    if(display_it)
    subplot(1,2,2); hold off;
    DisplayPoints2(Model,Scene);
    set(gca,'FontSize',16);
    title(sprintf('KC Value: %f',ER(k+1)));
    drawnow;
    end
    
    if k>=k_flag+6
        if abs(ER(k+1)-ER(k-5))*10000000<=1 && abs(ER(k+1)-ER(k-8))*10000000<=1 && abs(ER(k+1)-ER(k-10))*10000000<=1 && k>200
            stop = 1;
        end
    end
end
% k
er=ER(1:k+1);
if not(arg.ReturnAll)
    TR = TR(:,:,k+1);
    TT = TT(:,:,k+1);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [match mindist] = match_bruteForce(q, p)
    m = size(p,2);
    n = size(q,2);    
    match = zeros(1,m);
    mindist = zeros(1,m);
    for ki=1:m
        d=zeros(1,n);
        for ti=1:3
            d=d+(q(ti,:)-p(ti,ki)).^2;
        end
        [mindist(ki),match(ki)]=min(d);
    end
    
    mindist = sqrt(mindist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [match mindist] = match_Delaunay(q, p, DT)
	match = transpose(nearestNeighbor(DT, transpose(p)));
	mindist = sqrt(sum((p-q(:,match)).^2,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [match mindist] = match_kDtree(~, p, kdOBJ)
	[match mindist] = knnsearch(kdOBJ,transpose(p));
    match = transpose(match);
function [match mindist] = match_kDtree2(q, p, params)
% 	[match mindist] = knnsearch(kdOBJ,transpose(p));
    [match,mindist] = flann_search(double(q),double(p),1,params);
%     match = transpose(match);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R,T] = eq_point(q,p,weights)

m = size(p,2);
n = size(q,2);

% normalize weights
weights = weights ./ sum(weights);

% find data centroid and deviations from centroid
q_bar = q * transpose(weights);
q_mark = q - repmat(q_bar, 1, n);
% Apply weights
q_mark = q_mark .* repmat(weights, 3, 1);

% find data centroid and deviations from centroid
p_bar = p * transpose(weights);
p_mark = p - repmat(p_bar, 1, m);
% Apply weights
%p_mark = p_mark .* repmat(weights, 3, 1);

N = p_mark*transpose(q_mark); % taking points of q in matched order

[U,~,V] = svd(N); % singular value decomposition

R = V*diag([1 1 det(U*V')])*transpose(U);

T = q_bar - R*p_bar;

function [T] = eq_point_T(q,p,weights)
% normalize weights
weights = weights ./ sum(weights);

% find data centroid and deviations from centroid
q_bar = q * transpose(weights);
% find data centroid and deviations from centroid
p_bar = p * transpose(weights);

T = q_bar - p_bar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R,T] = eq_plane(q,p,n,weights)

n = n .* repmat(weights,3,1);

c = cross(p,n);

cn = vertcat(c,n);

C = cn*transpose(cn);

b = - [sum(sum((p-q).*repmat(cn(1,:),3,1).*n));
       sum(sum((p-q).*repmat(cn(2,:),3,1).*n));
       sum(sum((p-q).*repmat(cn(3,:),3,1).*n));
       sum(sum((p-q).*repmat(cn(4,:),3,1).*n));
       sum(sum((p-q).*repmat(cn(5,:),3,1).*n));
       sum(sum((p-q).*repmat(cn(6,:),3,1).*n))];
   
X = C\b;

cx = cos(X(1)); cy = cos(X(2)); cz = cos(X(3)); 
sx = sin(X(1)); sy = sin(X(2)); sz = sin(X(3)); 

R = [cy*cz cz*sx*sy-cx*sz cx*cz*sy+sx*sz;
     cy*sz cx*cz+sx*sy*sz cx*sy*sz-cz*sx;
     -sy cy*sx cx*cy];
    
T = X(4:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R,T] = eq_lmaPoint(q,p)

Rx = @(a)[1     0       0;
          0     cos(a)  -sin(a);
          0     sin(a)  cos(a)];
      
Ry = @(b)[cos(b)    0   sin(b);
          0         1   0;
          -sin(b)   0   cos(b)];
      
Rz = @(g)[cos(g)    -sin(g) 0;
          sin(g)    cos(g)  0;
          0         0       1];

Rot = @(x)Rx(x(1))*Ry(x(2))*Rz(x(3));

myfun = @(x,xdata)Rot(x(1:3))*xdata+repmat(x(4:6),1,length(xdata));


options = optimset('Algorithm', 'levenberg-marquardt');
x = lsqcurvefit(myfun, zeros(6,1), double(p), double(q), [], [], options);


R = Rot(x(1:3));
T = x(4:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the RMS error between two point equally sized point clouds with
% point correspondance.
% ER = rms_error(p1,p2) where p1 and p2 are 3xn matrices.

function ER = rms_error(p1,p2)
dsq = sum(power(p1 - p2, 2),1);
ER = sqrt(mean(dsq));


function ER = rms_error2(p1,p2,a,weights)
if a==2
    ER=1;
%     dsq = sum(power(p1 - p2, 2),1);
%     ER=1/2*(mean(dsq))^2;
elseif a==0
    dsq = sum(power(p1 - p2, 2),1);
%     ER=log(1/2*(mean(dsq))^2+1);
    ER=2/(mean(dsq)+2);
elseif a==-inf
    dsq = sum(power(p1 - p2, 2),1);
% %     ER=1-exp(-((mean(dsq))^2/2));
    ER=exp(-(mean(dsq)/2));
%     size(weights)
%     size(p1')
%     MKDE = ComputeKDE22(p1',p2',weights);
%     ER = sum(MKDE);
else
    dsq = sum(power(p1 - p2, 2),1);
%     ER=abs(a-2)/a*((mean(dsq)/abs(a-2)+1)^(a/2)-1);
    ER=(mean(dsq)/abs(a-2)+1)^(a/2-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Converts (orthogonal) rotation matrices R to (unit) quaternion
% representations
% 
% Input: A 3x3xn matrix of rotation matrices
% Output: A 4xn matrix of n corresponding quaternions
%
% http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion

function quaternion = rmat2quat(R)

Qxx = R(1,1,:);
Qxy = R(1,2,:);
Qxz = R(1,3,:);
Qyx = R(2,1,:);
Qyy = R(2,2,:);
Qyz = R(2,3,:);
Qzx = R(3,1,:);
Qzy = R(3,2,:);
Qzz = R(3,3,:);

w = 0.5 * sqrt(1+Qxx+Qyy+Qzz);
x = 0.5 * sign(Qzy-Qyz) .* sqrt(1+Qxx-Qyy-Qzz);
y = 0.5 * sign(Qxz-Qzx) .* sqrt(1-Qxx+Qyy-Qzz);
z = 0.5 * sign(Qyx-Qxy) .* sqrt(1-Qxx-Qyy+Qzz);

quaternion = reshape([w;x;y;z],4,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Converts (unit) quaternion representations to (orthogonal) rotation matrices R
% 
% Input: A 4xn matrix of n quaternions
% Output: A 3x3xn matrix of corresponding rotation matrices
%
% http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#From_a_quaternion_to_an_orthogonal_matrix

function R = quat2rmat(quaternion)
q0(1,1,:) = quaternion(1,:);
qx(1,1,:) = quaternion(2,:);
qy(1,1,:) = quaternion(3,:);
qz(1,1,:) = quaternion(4,:);

R = [q0.^2+qx.^2-qy.^2-qz.^2 2*qx.*qy-2*q0.*qz 2*qx.*qz+2*q0.*qy;
     2*qx.*qy+2*q0.*qz q0.^2-qx.^2+qy.^2-qz.^2 2*qy.*qz-2*q0.*qx;
     2*qx.*qz-2*q0.*qy 2*qy.*qz+2*q0.*qx q0.^2-qx.^2-qy.^2+qz.^2];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Least squares normal estimation from point clouds using PCA
%
% H. Hoppe, T. DeRose, T. Duchamp, J. McDonald, and W. Stuetzle. 
% Surface reconstruction from unorganized points. 
% In Proceedings of ACM Siggraph, pages 71:78, 1992.
%
% p should be a matrix containing the horizontally concatenated column
% vectors with points. k is a scalar indicating how many neighbors the
% normal estimation is based upon.
%
% Note that for large point sets, the function performs significantly
% faster if Statistics Toolbox >= v. 7.3 is installed.
%
% Jakob Wilm 2010

function n = lsqnormest(p, k)
m = size(p,2);
n = zeros(3,m);

v = ver('stats');
if str2double(v.Version) >= 7.5 
    neighbors = transpose(knnsearch(transpose(p), transpose(p), 'k', k+1));
else
    neighbors = k_nearest_neighbors(p, p, k+1);
end

for i = 1:m
    x = p(:,neighbors(2:end, i));
    p_bar = 1/k * sum(x,2);
    
    P = (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k)); %spd matrix P
    %P = 2*cov(x);
    
    [V,D] = eig(P);
    
    [~, idx] = min(diag(D)); % choses the smallest eigenvalue
    
    n(:,i) = V(:,idx);   % returns the corresponding eigenvector    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program to find the k - nearest neighbors (kNN) within a set of points. 
% Distance metric used: Euclidean distance
%
% Note that this function makes repetitive use of min(), which seems to be
% more efficient than sort() for k < 30.

function [neighborIds neighborDistances] = k_nearest_neighbors(dataMatrix, queryMatrix, k)

numDataPoints = size(dataMatrix,2);
numQueryPoints = size(queryMatrix,2);

neighborIds = zeros(k,numQueryPoints);
neighborDistances = zeros(k,numQueryPoints);

D = size(dataMatrix, 1); %dimensionality of points

for i=1:numQueryPoints
    d=zeros(1,numDataPoints);
    for t=1:D % this is to avoid slow repmat()
        d=d+(dataMatrix(t,:)-queryMatrix(t,i)).^2;
    end
    for j=1:k
        [s,t] = min(d);
        neighborIds(j,i)=t;
        neighborDistances(j,i)=sqrt(s);
        d(t) = NaN; % remove found number from d
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Boundary point determination. Given a set of 3D points and a
% corresponding triangle representation, returns those point indices that
% define the border/edge of the surface.

function bound = find_bound(pts, poly)

%Correcting polygon indices and converting datatype 
poly = double(poly);
pts = double(pts);

%Calculating freeboundary points:
TR = TriRep(poly, pts(1,:)', pts(2,:)', pts(3,:)');
FF = freeBoundary(TR);

%Output
bound = FF(:,1);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 