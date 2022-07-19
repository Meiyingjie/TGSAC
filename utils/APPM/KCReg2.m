% function [param] = KCReg(M, S, h, display,motion);
%M is a N by 2 array containing 2D model points
%S is a M by 2 array containing 2D Scene points
%h is the "bandwith"
%display: display the intermediate steps or not. default is not display
%motion:  the transformation model, can be
%         euclidean
%         affine
%         projective;
%         Default motion model is euclidean;
function [Ricp, Ticp, ER] = KCReg2(M,S,config, h, display,flag,motion)

if nargin<3
    disp('Not enough input parameters');
    return;
end

if nargin<6
    display = 0;
    motion = 'euclidean';
end;

if nargin<7
    motion = 'euclidean';
end;

%set the global parameters that can be seen by "ComputeKC".
global display_it;
display_it = display;
global Scene;
Scene = S;
global Model;
Model = M;
global resolution;
resolution = h;
global min_val;
min_val = min(S);
global max_val;
max_val = max(S);
global SceneKDE;
% SceneKDE = ComputeKDE2(S);
SceneKDE =1;

if(display_it)
    figure(2);
    subplot(1,2,1); hold off;
    DisplayPoints2(Model,Scene);
    set(gca,'FontSize',16);
    title('Initial setup');
    drawnow;
end;


switch lower(motion)
    case 'euclidean'
%         opt = optimset('MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-6,'TolX',1e-10);
%         param = fminsearch('ComputeKC2',[0,0,0,0,0,0]',opt);
          if flag==1
%               [Ricp, Ticp, ER] = ICP_kdtree(M', S',config, 'Matching', 'kDtree','translation',true, 'WorstRejection', 0.2, 'iter', 1);
              [Ricp, Ticp, ER] = ICP_kdtree(M', S',config, 'Matching', 'kDtree', 'WorstRejection', 0.8, 'iter', 100);
          else
              [Ricp, Ticp, ER] = ICP_kdtree(M', S',config, 'Matching', 'kDtree', 'WorstRejection', 0.8, 'iter', 500);
          end
%           param=[Ricp, Ticp]
%           TMICP = Ricp*M' + repmat(Ticp, 1, size(M', 2));
% %           figure;
%           subplot(1,2,1); hold off;
%           DisplayPoints2(TMICP',M);
%           set(gca,'FontSize',16);
%           title('Initial setup');
%           drawnow;
    case 'affine'
        opt = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',1e-6,'TolX',1e-10);
        param = fminsearch('ComputeKC2',[1 0 0 1 0 0]',opt);        
    case 'projective'
        opt = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',1e-6,'TolX',1e-10);
        param = fminsearch('ComputeKC',[1 0 0 0 1 0 0 0 ]',opt);
        resolution = resolution/3;
        SceneKDE = ComputeKDE2(S);
        param = fminsearch('ComputeKC',param,opt);
    otherwise
        disp('Unknown motion type');
        return;
end;
