function [maps, mapeq, maps_b, mapeq_b] = genSDPConstraintMap2(nn,k, W)

    maps =  {}; %cell(2*n,1);     
    mapeq = {};
    eqCount = 1; ieqCount = 1;
    n=length(W);
%     Av(1:n) = 1;
%     size(Av)
%     for i=1:1
%         A = zeros(n,n);
%         A(i,:) = 1;
        

    M = sparse(n, n);
    M(n, :) = 1;
    maps{ieqCount} = M; ieqCount = ieqCount + 1;
%     M = sparse(n, n);
%     M(:, n) = 1;
%     maps{ieqCount} = M; 
    maps_b =1;
%     end
    
%     for i=1:n
%         A = zeros(n,n);
%         A(:,i) = 1;
%         Av = [A(:); 0];
%         
%         M = sparse(n*n+1, n*n+1);
%         M(n*n+1, :) = Av;
%         maps{ieqCount} = M; ieqCount = ieqCount + 1;
%     end
%     maps_b = ones (2*nn, 1);
     
 %----------------------------------------------------
    
    % 1^TX1 = k
    M = sparse(n, n);
    M(n,1:(n-1)) = 1;
    mapeq{eqCount} = M; 
    mapeq_b = k;

end