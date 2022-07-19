clc;clear;close all
base=importdata('F:\新桌面\资料1\数字孪生\全新文献\点云\程序\数据\facade\groundtruth\facade.txt');
tran=importdata('C:\Users\23560\Desktop\pointset\数据\面积facade点云坐标.txt');
num=length(base)/4;
for i=1:num
      flag1=[tran(3+(i-1)*4,1),tran(3+(i-1)*4,2),tran(1+(i-1)*4,1)];
      flag2=[base(3+(i-1)*4,1),base(3+(i-1)*4,2),base(1+(i-1)*4,1)];
      R1_record(:,i)=solve_R(flag1);
      R2_record(:,i)=solve_R(flag2);
      T1_record(:,i)=tran(1+(i-1)*4:3+(i-1)*4,4);
      T2_record(:,i)=base(1+(i-1)*4:3+(i-1)*4,4);
end
for j1=1:3
    for j2=1:num
        R_erro(j1,j2)=abs(R1_record(j1,j2)-R2_record(j1,j2));
        T_erro(j1,j2)=abs(T1_record(j1,j2)-T2_record(j1,j2));
    end
end
sum(sum(R_erro*180/pi))/30
sum(sum(T_erro))/30
% for i=1:16
%     for j=1:4
%         erro(i,j)=abs(tran(i,j)-base(i,j));
%     end
% end
% sum(sum(erro))/64
% sin(30*pi/180)

function pp=solve_R(flag)
t=[0;0;0]*pi/180;
pp=fsolve(@(x)myfun(x,flag),t);
end
% erro=abs(tran-base);
% sum(erro)/64
function F = myfun(x,flag)
F = [-sin(x(2))-flag(1);
      cos(x(2))*sin(x(1))-flag(2);
      cos(x(3))*cos(x(2))-flag(3)];
end