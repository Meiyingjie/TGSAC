% a=randn(5,6);
% [m1,idx]=max(a);
% % idy=1:size(a,2);
% [m2,idx2]=max(m1);
% a(idx(idx2),idx2)
% [m3,idx3]=max(max(a))

clc;clear;close all
gai=0.99;
lv=0.01;

% lv=3/64;
% gai=1-lv;
uyuy=0;
m=64;
n=3;
% nchoosek(pp,3)
for i=1:n
    yuyu=vpa(gai^(m-i+1)*lv^(i-1)*nchoosek(m,i-1),10);
    
    uyuy=yuyu+uyuy;
%     yuyu=vpa(gai^(n-i+1)*lv^(i-1),10);
%     uyuy=yuyu+uyuy;
end
uyuy2=lv^32;
T1=log(0.01)/log(uyuy)
T2=log(0.01)/log(1-uyuy2)
lv*(gai^3+gai^2*lv+lv^2*gai)+gai^4



%%
% q=[0.00449209 0.38422 -0.00976512 0.923179];
% t=[-0.00646017 -1.36122e-05 -0.0129064];
% % R=[1-2*q(3)^2-2*q(4)^2 2*q(2)*q(3)+2*q(1)*q(4) 2*q(2)*q(4)-2*q(1)*q(3)
% %    2*q(2)*q(3)-2*q(1)*q(4) 1-2*q(2)^2-2*q(4)^2 2*q(3)*q(4)+2*q(1)*q(2)
% %    2*q(2)*q(4)+2*q(1)*q(3) 2*q(3)*q(4)-2*q(1)*q(2) 1-2*q(2)^2-2*q(3)^2];
% R=[2*q(3)^2+2*q(4)^2-1 -2*q(2)*q(3)-2*q(1)*q(4) -2*q(2)*q(4)+2*q(1)*q(3)
%    -2*q(2)*q(3)+2*q(1)*q(4) 2*q(2)^2+2*q(4)^2-1 -2*q(3)*q(4)-2*q(1)*q(2)
%    2*q(2)*q(4)+2*q(1)*q(3) 2*q(3)*q(4)-2*q(1)*q(2) 1-2*q(2)^2-2*q(3)^2];
% trans=[R [t(1);t(2);t(3)]]


