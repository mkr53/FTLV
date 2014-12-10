
function [W,H,phi] = ACLS()%V,C, sunsky)
    %Here we use alternating constrained least squares to decompose the nxm
    %matrix V into a nxk matrix W and a kxm matrix H. 
    %The concept is to minimize the euclidean error: 
    %||V-WH||^2 = sigma(Vij - (WH)ij)2
    
    % set up paths to VLFeat functions. 
    % See http://www.vlfeat.org/matlab/matlab.html for VLFeat Matlab documentation
    % This should work on 32 and 64 bit versions of Windows, MacOS, and Linux
        run('/Users/elisecotton/vlfeat-0.9.19/toolbox/vl_setup')
        
    %small matrices for testing 
   
    fun = @(x)x*x;
    x = 100;
    z = 50;
    temp = zeros(x,z);
    begin_matrix = linspace(1,10,x);
    for g = 1:z
        temp(:,g) = begin_matrix;
    end
    
    %V = arrayfun(fun, temp);
    V = randi( [0 256], 10, 50);
    [n,m] = size(V);
    C = randi([0 1], n, m);
    sunsky= 'sun'
    
    
    k = 1; %k is between 1 and 5. paper uses 1
    
    fprintf('initializing W and H \n')
    [W,H] = initializeACLS(V,k);  
    
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
    %test against other decomposer    
    NewA = V .* C;
    [Wcoeff,Hbasis,numIter,tElapsed,finalResidual]=wnmfrule(NewA,1);
    
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
    myfun = @(ph,t)H(round(t + ph));
    
    %phi = calculatePhi(phi,C,V,H,n,m);
     phi = calculatePhi(V,H,n);


    
    fprintf('solving for W and H \n');
    %iterations
    its = 1;
    
    %we will alternate between linear squares solves of W and H:
    for j = 1:its
        fprintf('optimizing W \n');
        %for each row of W, solve LCLS: M = H' d = v' x = w'        
        for i = 1:n
           %shift entire matrix H by corresponding phi value
           %H_shift = circshift(H,[phi(i,:),0]);
           H_shift = shiftByPhi(H, phi(i,:),m);
           %to include confidence, M = C_i * H_shift' and d = C_i * V_i
           M = double(C(i,:)' .* (H_shift'));
           d = double((C(i,:) .* V(i,:))');
           
           x = lsqnonneg(M,d);        
           W(i,:) = x';
        end
        
        fprintf('optimizing H \n');
        %for each column of H, solve LCLS: M = W d = v x = h
        for i = 1:m
            %shift entire column of V by -phi 
            V_shift = zeros(n,1);
            for d = 1:n 
                %new column to grab data from(clamped between 1 and m) 
                k = max(min((i + phi(d)),m),1);
                V_shift(d) = V(d,k);
            end
           
            %to include confidence, M = C_i * W and d = C_i * V_i
            M = double(C(:,i) .* W);
            d = double(C(:,i) .* V_shift);
            
            x = lsqnonneg(M,d);
            H(:,i) = x;
        end
        
        fprintf('optimizing phi \n');
        %now we solve for phi
         phi = calculatePhi(V,H,n);
        
    end
    
    i = round(n / 2);
    figure
    frames = 1:m;
    colors = ['c','y','g','b'];
    %shifted_frames = frames + phi(i,:);
        hold on;
        
    for it = 1:2
        shifted_frames = frames - phi(i,:);

        phi(i,:)
        plot(frames, C(i,:) .* V(i,:),colors(2 + it));
        
        plot(shifted_frames, C(i,:) .* V(i,:),colors(it));
        i = i + round(n/8);
    end
            plot(frames, H,'r');

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

function phi_n = calculatePhi(V,H,n)
phi_n = zeros(n,1);
    for i = 1:n
        [r,lag] = xcorr(V(i,:), H);
        [~,I] = max(abs(r));
        lagDiff = lag(I);
        phi_n(i) = -1 * lagDiff;
    end
end

function phi_n = calculatePhi_o(phi,C,V,H,n,m)
    %used for phi estimation
    x_data = 1:m;
    my_fun = @(ph,t)H(round(t + ph));
    
    for i = 1:n
        options = optimoptions('lsqcurvefit', 'Display', 'off');
        ph = lsqcurvefit(my_fun, 0, x_data, (C(i,:) .* V(i,:)),[],[],options);
        phi(i) = ph;
        
    end
    
    phi_n = phi;
end

function H_shift = shiftByPhi(H,phi,f)
H_shift = zeros(1,f);
shifted_indices = max(1 + phi , 1) : min(f + phi, f);
H_shift(shifted_indices) = H(:,1:numel(shifted_indices));
end