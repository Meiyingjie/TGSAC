%compute the Kernel Density Estimation defined by the N by 2 array P
function KDE = ComputeKDE22(P,PP,weights)

% diff=P-PP;
% size(-diff.*(diff))
% % size(weights./(sum(weights)))
% KDE=weights/(sum(weights))*exp(-diff.*diff)/2;

for flag=1:min(size(P,1),size(PP,1))
    diff=P(flag,:)-PP(flag,:);
%     KDE(flag)=exp(-diff*diff'/2)*norm(diff);
    KDE(flag)=exp(-diff*diff'/2)*norm(diff);
end
% size(KDE')
% KDE=weights/(sum(weights))*KDE';
