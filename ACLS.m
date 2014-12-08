
function [W,H,phi] = ACLS(V,C, sunsky)
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
    

    fprintf('initializing W and H \n')
    [W,H] = initializeACLS(V,k);  
    
    
    
    %not necessary in lsqnonneg
    %A = -1 * eye(k);
    %b = zeros(k,1);
    
    if strcmp(sunsky,'sun')
        [W,H,phi] = sun_solve(V,C,W,H);
    elseif strcmp(sunsky,'sky')
        [W,H,phi] = sky_solve(V,C,W,H);
    else
        ME = MException('Invalid input for sunsky',sunsky);
        throw(ME)
    end
end
function [W_f,H_f,phi] = sky_solve(V,C,W,H)
    %n = # pixels ;  m = # frames 
    [n,m] = size(V);
    
    fprintf('solving for W and H \n');
    %iterations
    its = 1;
    
    %we will alternate between linear squares solves of W and H:
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
    
    W_f = W;
    H_f = H;
    phi = zeros(n,1);
end

function [W_f,H_f,phi] = sun_solve(V,C,W,H)
    %n = # pixels ;  m = # frames 
    [n,m] = size(V);
    %initialize phi, our shift map if sunsky = 'sun'
    phi = zeros(n,1);
    %used for phi estimation
    x_data = 1:m;
    myfun = @(ph,t)H(t) + ph;
    

    fprintf('solving for W and H \n');
    %iterations
    its = 1;
    
    %we will alternate between linear squares solves of W and H:
    for j = 1:its
        %for each row of W, solve LCLS: M = H' d = v' x = w'        
        for i = 1:n
           %shift entire matrix H by corresponding phi value
           H_shift = H + phi(i,:);
           %to include confidence, M = C_i * H_shift' and d = C_i * V_i
           M = double(C(i,:)' .* (H_shift'));
           d = double((C(i,:) .* V(i,:))');
           
           x = lsqnonneg(M,d);        
           W(i,:) = x';
        end

        %for each column of H, solve LCLS: M = W d = v x = h
        for i = 1:m
            %shift entire column of V by -phi 
            V_shift = V(:,i) - phi;
            %to include confidence, M = C_i * W and d = C_i * V_i
            M = double(C(:,i) .* W);
            d = double(C(:,i) .* V_shift);
            
            x = lsqnonneg(M,d);
            H(:,i) = x;
        end
        %now we solve for phi
        %for i = 1:n
            i = 50000;
            frames = 1:m; 
            
            figure
            subplot(1,1,1)
            hold on;
            plot(frames, (C(i,:) .* V(i,:)),'m');
            plot(frames, H,'r');

            hold off;
            
            options = optimoptions('lsqcurvefit', 'Display', 'off');
            ph = lsqcurvefit(myfun, 0, x_data, (C(i,:) .* V(i,:)),[],[],options);
            phi(i) = ph;
            ph
        %end
    end
    
    W_f = W;
    H_f = H;
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