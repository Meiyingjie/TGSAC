% given the transformation parameter param, compute the KC value;
function KCVal = ComputeKC2(R,T,weights)
global Scene;
global Model;
global SceneKDE;
global display_it;
global resolution;

PT = TransformPoint2(R,T,Model);
MKDE = ComputeKDE22(Scene,PT,weights);
% MKDE.*SceneKDE
KCVal = -sum(MKDE);
% MKDE = ComputeKDE2(PT);
% % MKDE.*SceneKDE
% KCVal = -sum(sum(MKDE.*SceneKDE));
% sum(sum(MKDE.*SceneKDE))

%The following for display purpose only.
if(display_it)
    subplot(1,2,2); hold off;
    DisplayPoints2(PT,Scene);
    set(gca,'FontSize',16);
    title(sprintf('KC Value: %f',-KCVal));
    drawnow;
end
