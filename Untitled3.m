base=[0.1125   -0.9588    0.2609 -0.12392735
    0.8369    0.2330    0.4953 -0.936076188
   -0.5357    0.1626    0.8286 0.52748802];

tran=[[0.0487590751239061,-0.939867601676410,0.338040595065332,-0.199614660100554;0.746178220852821,0.259256390305821,0.613191802628287,-1.06062269944190;-0.663958293357369,0.222339864632177,0.713949836667358,0.829726706819115]];

R_base=base(1:3,1:3);
R_tran=tran(1:3,1:3);
T_base=base(1:3,4);
T_tran=tran(1:3,4);
%       flag=(trace(R_base*R_tran')-1)/2;
%       erro_R(i)=solve_R(flag);     
erro_R=acosd((trace(R_base*R_tran')-1)/2);
erro_T=norm(T_base-T_tran);

erro_R
erro_T
function pp=solve_R(flag)
t=0.2*pi/180;
pp=fsolve(@(x)myfun(x,flag),t);
end
% erro=abs(tran-base);
% sum(erro)/64
function F = myfun(x,flag)
F = cos(x)-flag;
end