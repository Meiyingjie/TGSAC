% Given two point set M and B, find the set of correspondences using SDP

function [R,t, corrM, corrB] = sdpReg2(M, B, config)
    n=length(M);
    % Generate Weight Matrix
    corrM = []; corrB = []; R=[]; t = [];    
%     [W, acceptedPairs] = generateWeightMatrix(M, B, config); 
    [W, acceptedPairs] = generateWeightMatrix3(M, B, config);
%     acceptedPairs
%     W
%     size(W)
%     acceptedPairs=length(index);
%     size(Tx)
%     acceptedPairs
    if (acceptedPairs > config.nareaPairThresh)
        
        k = config.k;
        
%         [X] = solveSDP(W,k,n); 
%         [X] = solveSDP2(W,k,n);
        [X] = solveSDP3(W,config);
        [sM, sB, corrB] = getCorrespondences2(X,M,B);
   
        [R, t] = Kabsch(sM, sB);
        
    end
    
    

end

