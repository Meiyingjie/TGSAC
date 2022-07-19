% Main SDRSAC Algorithm
% Input: 
%       M, B: Input point clouds, each is a 3xN matrix  
%       config: Contains all parameters required to run the algoirthm.
%               See readConfig.m for more info    
% Output: a structure out containing output variables:
%         out.bestR: Best rotation matrix
%         out.bestT: Best translation vector


function [out] = SDRSAC2(M, B, MM, BB, config)  
    add_dependencies;   
    [Ricp1, Ticp1] =KCReg2(BB', MM',config, 2,0,0);
%     Ricp1=[1 0 0;0 1 0;0 0 1 ];
%     Ticp1=[0;0;0];
    TMICP1 = Ricp1*M + repmat(Ticp1, 1, size(M, 2)); 
%     TMICP2 = Ricp1*MM + repmat(Ticp1, 1, size(MM, 2)); 
    ps = 0.99;
    M=TMICP1;
    ceshi=TMICP1';
%     MM=TMICP2;
    iter=0;
    T_max = 1e10;
%     TR=ones(3,4);
    figure;
    plotPointClouds(M, B, 'b.','r.')
    % Using KDTree for quick computing of consensus size
    B_Tree = KDTreeSearcher(B');
    inls_icp = countCorrespondences(M, B_Tree, config.epsilon);
%     seg_size=max(floor(max(length(M),length(B))/200),2);
    % Prepare sampling
    n = config.pointPerSample;  
    maxInls = inls_icp;
    fprintf('≥ı ºªØ consensus size: %d\n', maxInls);%
    fprintf('--------------');
    bestR = {};
    bestT = {};
%     lap=[];
    bestR{1}=Ricp1;
    bestT{1}=Ticp1;
    ER=0;
    pI =inls_icp./size(M,2);
    lap(1)=20;
    B_to_sample = B;
    % Start the sampling iterations
    stop = false;
    stop_flag=0;
    record=1;
%     TMICP=B;
    while (iter < config.maxIter && ~stop)        
        if(length(B_to_sample)<=length(B)*0.95)
            B_to_sample=B;
%         elseif length(B_to_sample)<=length(B)*0.99
%             B_to_sample=TMICP;
        end  
%         if(length(B_to_sample)<=length(B)*0.8)
%             B_to_sample=TMICP;
%         end
        tic;
        [id1,id2]=suoyin(M,B_to_sample,config,3);
%          B_to_sample = B;
%         idxM = randsample(size(M,2), n);
%         m = M(:, idxM);     
        idxM=id1;
        m = M(:, unique(idxM)); 
%         scount = 0;
%         while (size(B_to_sample,2) > n*2 && scount < 1)
%             scount = scount + 1;
            fprintf('Current B size: %d\n', size(B_to_sample,2));
%             BB=B_to_sample;
%             [~,id2]=suoyin(M,B_to_sample,config,5);
%             idxB = randsample(size(B_to_sample,2), n);
            idxB=id2;
            b = B_to_sample(:, unique(idxB));  
%             size(id2)
%             b = B_to_sample(:, id2); 
            B_to_sample(:,unique(idxB)) = [];
            
            % Solve SDP
%             [Rs,ts, ~, corrB] = sdpReg(m, b, config);
          
            [Rs,ts, ~, corrB] = sdpReg2(m, b, config);
            toc;
%             corrB
%             figure;
%             plotPointClouds(M, B_to_sample, 'b.','r.')
%             hold on;
%             plotPointClouds(m, b, '+w','+k')
            if ~isempty(corrB)
                TM = Rs*M + repmat(ts, 1, size(M,2));
               
                % Conducting ICP
%                 [Ricp, Ticp] = icp(B, TM, 'Matching', 'kDtree', 'WorstRejection', 0.2, 'iter', 500);  
                [Ricp, Ticp, ER] =KCReg2(B', TM',config, 2,0,0);
%                   [Ricp, Ticp] = TrICP(B, TM, 100, 1);
                TMICP = Ricp*Rs*M + repmat(Ricp*ts + Ticp, 1, size(M, 2));   
%                 figure;
%                 set(gca,'FontSize',12);
%                 plot3(B(1,:),B(2,:),B(3,:),'k+');
%                 hold on;
%                 plot3(TMICP(1,:),TMICP(2,:),TMICP(3,:),'go');
%                 axis equal;
                inls_icp = countCorrespondences(TMICP, B_Tree, config.epsilon);   
%                 B_to_sample=TMICP;
                if inls_icp > maxInls
                    stop_flag=0;
                    record=record+1;
                    maxInls = inls_icp;
%                     bestR{record} = Ricp*Rs;
%                     bestT{record} = Ricp*ts + Ticp;  
                    bestR{record}=Ricp*Rs*bestR{1};
                    bestT{record}=Ricp*Rs*bestT{1} + Ricp*ts + Ticp;  
                    
                    fprintf('Best-so-far consensus size: %d\n', maxInls);%
                    fprintf('--------------');
%                     iter_flag=1;
%                     if(maxInls>=(length(B)*0.2))
%                     B_to_sample=TMICP;
%                     end
                    % For debugging purpose:
                     %close all; plotPointClouds(B, TMICP, 'b.','r.');
                    
                    % Compute stopping criterion
                    pI = maxInls./size(M,2);
                    T_max = log(1-ps)./log(1-pI^config.k);
                else
                    stop_flag=stop_flag+1;
                    stop_flag
                    if(stop_flag==20 && config.overlap==0)
                        stop = true;
                    end
                end                
            end 
            uyuy=0;
            for i=1:3
                yuyu=vpa((1-pI)^(64-i+1)*pI^(i-1)*nchoosek(64,i-1),10);
                uyuy=yuyu+uyuy;
            end
            if uyuy==1
                uyuy=0.99999999;
            end
%             uyuy
%             uyuy*pI
            iter_max=log(0.01)/log(1-(1-uyuy)*(pI));
            iter = iter+1;
%             iter_flag=iter_flag+1;
%             abs(ER(end)+0.5)*10000000
%             oo2=ER(end)
%             roundn(ER(end)*10000,1)
%             (roundn(ER(end)*10000,1)==-5000)
            % Stopping criterion is satisfied (minimum of 5 iterations)
            if T_max==-Inf
                T_max=1000;
            end
            if iter_max==-Inf
                iter_max=500;
            end
            if config.overlap~=0
                lap(record)=abs(pI-config.overlap);
                [~,uu]=min(lap);
                if abs(pI-config.overlap)<=0.1*config.overlap || abs(ER(end)+1.0)*10000000<=1 ||  stop_flag>5&&iter>1
%                     fprintf(lap,'/n');
%                     fprintf(uu,'/n');
                    bestR{end}=bestR{uu};
                    bestT{end}=bestT{uu};
                    stop = true;
                end
%             elseif abs(ER(end)+0.5)*10000000<=1 || (iter >= T_max && iter >=5 )||  stop_flag>5
            elseif abs(ER(end)+1.0)*10000000<=1 || (stop_flag >= iter_max) ||pI>=0.5||iter>500&&iter>1
%             if iter >= 100 && iter >=5 
                 stop = true;
            end
%         end
        
    end
    if length(bestR)>1
        tran=bestR{end}*MM + bestT{end};
%         tran=bestR{end}*Ricp1*MM + bestR{end}*repmat(Ticp1, 1, size(MM, 2)) + bestT{end};
%         best_R=bestR{end}*Ricp1;
%         best_T=bestR{end}*Ticp1 + bestT{end};
    else
        tran=Ricp1*MM + repmat(Ticp1, 1, size(MM, 2));
    end
%     out.ER=ER;
    out.inls = maxInls;    
    out.iter = iter;    
    out.R = bestR;
    out.T = bestT;
    out.best_RT= [bestR{end},bestT{end}];
    out.ori=BB';
    out.TMICP=tran';
    figure;
%     plotPointClouds(BB, tran, 'b.','r.')
    plotPointClouds(BB, tran)
%     figure;
%     pcshowpair(B, bestR*M + bestT);
end 
    