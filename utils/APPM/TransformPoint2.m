%2D transform of point set
% param = (dx,dy,theta);
% transformation is done by first shift the point center back origin,
% rotate, then shift back
function PT = TransformPoint2(R,T,P)
PT = (R* P')';
% r=rt1*rt2*rt3
PT(:,1) = PT(:,1) + T(1);
PT(:,2) = PT(:,2) + T(2);
PT(:,3) = PT(:,3) + T(3);
% switch length(param)
%     case 6
%         param
%         % Euclidean motion
%         center = mean(P);
%         rt1 = [1 0 0; 0, cos(param(4)), -sin(param(4)); 0 sin(param(4)), cos(param(4))];
%         rt2 = [cos(param(5)), 0, sin(param(5)); 0 1 0; -sin(param(5)),0, cos(param(5))];
%         rt3 = [cos(param(6)), -sin(param(6)), 0; sin(param(6)), cos(param(6)), 0; 0 0 1];
%         PT = (rt1*rt2*rt3 * [P(:,1)-center(1) P(:,2)-center(2) P(:,3)-center(3)]')';
%         r=rt1*rt2*rt3
%         PT(:,1) = PT(:,1) + param(1)+center(1);
%         PT(:,2) = PT(:,2) + param(2)+center(2);
%         PT(:,3) = PT(:,3) + param(3)+center(3);
%         center1 = mean(PT)
%     case 6
%         % Affine motion
%         M = reshape(param,2,3);
%         PT = (M * [P'; ones(1,size(P,1))])';
%     case 8
%         % Projective motion
%         M = reshape([param;1],3,3);
%         PT2 = (M * [P'; ones(1,size(P,1))])';
%         % make sure not divided by zero;
%         PT2(:,3) = sign(PT2(:,3)) .* max(abs(PT2(:,3)), 1e-10);
%         PT = zeros(size(P,1),2);
%         PT(:,1) = PT2(:,1)./PT2(:,3);
%         PT(:,2) = PT2(:,2)./PT2(:,3);
%     otherwise
%         disp('Unknown motion type');
%         error('','unknown motion');
end
