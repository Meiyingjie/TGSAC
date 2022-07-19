% Generate weight matrix for the task of quadratic matching
function  [W, acceptedPairs] = generateWeightMatrix3(M, B, config)
    spc1=MyCrust(double(M'));
    spc2=MyCrust(double(B'));
    N1=length(spc1);
    N2=length(spc2);
%     W = -10000*ones(max(N1,N2), max(N1,N2));
    W = -10000*ones(max(N1,N2), max(N1,N2));
    acceptedPairs = 0;
    % Loop through all possible pairs
    for t1=1:N1
        for t2=1:N2
            p=spc1(t1,:);
            q=spc2(t2,:);
            if max(p)>=length(M)
                [index1,index11,~]=find(p>=length(M));
                p(index1,index11)=length(M);
            end
            if max(q)>=length(B)
                [index2,index22,~]=find(q>=length(B));
                q(index2,index22)=length(B);
            end
            diff1=area(M,p(1),p(2),p(3),config);
            diff2=area(B,q(1),q(2),q(3),config);
            du1=angle(M,p(1),p(2),p(3),config);
            du2=angle(B,q(1),q(2),q(3),config);
            if diff1~=0&&diff2~=0
                dDiff = abs(diff1-diff2) ;
                du_cha1=abs(du1(1)-du2(1));
                du_cha2=abs(du1(2)-du2(2));
            end
            if dDiff <= config.dDiffThresh1&&du_cha1<=10&&du_cha2<=10
                W (t1,t2) = exp(-dDiff);  % can also use 1./dDiff; %
                acceptedPairs = acceptedPairs + 1; 
            end
%             if diff1~=0&&diff2~=0
%                 dDiff =  abs(diff1-diff2) ;
%             end
%             if dDiff <= config.dDiffThresh1
%                 W (t1,t2) = exp(-dDiff);  % can also use 1./dDiff; %
%                 acceptedPairs = acceptedPairs + 1; 
%             end
        end
    end
end
    
