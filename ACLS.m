
function [W,H] = ACLS()
    %Here we use alternating constrained least squares to decompose the nxm
    %matrix V into a nxk matrix W and a kxm matrix H. 
    %The concept is to minimize the euclidean error: 
    %||V-WH||^2 = sigma(Vij - (WH)ij)2
    
    % set up paths to VLFeat functions. 
    % See http://www.vlfeat.org/matlab/matlab.html for VLFeat Matlab documentation
    % This should work on 32 and 64 bit versions of Windows, MacOS, and Linux
        run('/Users/elisecotton/vlfeat-0.9.19/toolbox/vl_setup')
        
    load('imagematrix.mat', 'images');
    V = images;
    %V = randi( [0 256], 100, 500);
    k = 1; %k is between 1 and 5
    [n,m] = size(V);
    
    [W,H] = initializeACLS(V,k);
    
    %we will alternate between linear squares solves of W and H:
    %x = lsqlin(C,d,A,b)solves the linear system C*x = d in the 
    %least-squares sense, subject to A*x ? b.
    
    A = -1 * eye(k);
    b = zeros(k,1);
    
    its = 1;
    
    distances = zeros(its,1);
    
    for j = 1:its
        %for each row of W, solve LCLS: C=H' d = v' x = w'
        C = double(H');
        for i = 1:n
            d = double(V(i,:)');
            x = lsqlin(C,d,A,b);        
            W(i,:) = x';
        end

        %for each column of H, solve LCLS: C=W d=v x=h
        C = double(W);
        for i = 1:m
            d = double(V(:,i));
            x = lsqlin(C,d,A,b);
            H(:,i) = x;
        end

        %check distance
        WH = W*H;
        A2 = reshape(WH',numel(size(WH)),1); % makes column vectors
        B2 = reshape(V',numel(size(V)),1);

        dist = sqrt(dot(A2-B2,A2-B2));
        distances(j) = dist;
    end
    
    distances
end

function [W,H] = initializeACLS(V,k)
    %initialization:
    % first, run k-means clustering with k = 20. then, populate our W and H
    % matrices with random cluster centers 
    vocab_size = 20;
    
    [n,m] = size(V);
    
    [centers, assignments] = vl_kmeans(single(V), vocab_size);
    size(centers)
    
    centers_vector = reshape(centers,1,numel(centers));

    W = reshape(datasample(centers_vector,(n*k)), n, k);
    H = reshape(datasample(centers_vector, k*m),k,m);
end 