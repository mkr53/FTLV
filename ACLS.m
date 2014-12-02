
function [W,H] = ACLS(V,C, sunsky)
    %Here we use alternating constrained least squares to decompose the nxm
    %matrix V into a nxk matrix W and a kxm matrix H. 
    %The concept is to minimize the euclidean error: 
    %||V-WH||^2 = sigma(Vij - (WH)ij)2
    
    % set up paths to VLFeat functions. 
    % See http://www.vlfeat.org/matlab/matlab.html for VLFeat Matlab documentation
    % This should work on 32 and 64 bit versions of Windows, MacOS, and Linux
        run('/Users/elisecotton/vlfeat-0.9.19/toolbox/vl_setup')
        
    %small matrices for testing 
    %V = randi( [0 256], 100, 500);
    %C = randi([0 1], n, m);
    
    
    k = 1; %k is between 1 and 5. paper uses 1
    [n,m] = size(V);

    fprintf('initializing W and H \n')
    [W,H] = initializeACLS(V,k);  
    
    %we will alternate between linear squares solves of W and H:
    %x = lsqlin(C,d,A,b)solves the linear system C*x = d in the 
    %least-squares sense, subject to A*x ? b.
    
    %not necessary in lsqnonneg
    %A = -1 * eye(k);
    %b = zeros(k,1);
    
    fprintf('solving for W and H \n');
    %iterations
    its = 1;
    
    for j = 1:its
        %for each row of W, solve LCLS: M = H' d = v' x = w'        
        for i = 1:n
           %to include confidence, M = C_i * H' and d = C_i * V_i
           M = double(C(i,:)' .* (H'));
           d = double((C(i,:) .* V(i,:))');
           
           x = lsqnonneg(M,d);        
           W(i,:) = x';
        end

        %for each column of H, solve LCLS: M = W d = v x = h
        for i = 1:m
            %to include confidence, M = C_i * W and d = C_i * V_i
            M = double(C(:,i) .* W);
            d = double(C(:,i) .* V(:,i));
            
            x = lsqnonneg(M,d);
            H(:,i) = x;
        end
   
    end
     
end

function [W,H] = initializeACLS(V,k)
    %initialization:
    % first, run k-means clustering with k = 20. then, populate our W and H
    % matrices with random cluster centers 
    vocab_size = 20;
    
    [n,m] = size(V);
    
    [centers, assignments] = vl_kmeans(single(V), vocab_size);
    
    centers_vector = reshape(centers,1,numel(centers));

    W = double(reshape(datasample(centers_vector,(n*k)), n, k));
    H = double(reshape(datasample(centers_vector, k*m),k,m));
end 