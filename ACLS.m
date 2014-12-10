
function [W,H,phi,H_shifted] = ACLS(V,C, sunsky)
    %Here we use alternating constrained least squares to decompose the nxm
    %matrix V into a nxk matrix W and a kxm matrix H. 
    %The concept is to minimize the euclidean error: 
    %||V-WH||^2 = sigma(Vij - (WH)ij)2
    
    % set up paths to VLFeat functions. 
    % See http://www.vlfeat.org/matlab/matlab.html for VLFeat Matlab documentation
    % This should work on 32 and 64 bit versions of Windows, MacOS, and Linux
        run('/Users/elisecotton/vlfeat-0.9.19/toolbox/vl_setup')
    
        %{
    %small matrices for testing 
    V = randi( [0 256], 10, 50);
    [n,m] = size(V);
    C = randi([0 1], n, m);
    sunsky = 'sun'
    %}
    
    k = 1; %k is between 1 and 5. paper uses 1
    
    fprintf('initializing W and H \n')
    [W,H] = initializeACLS(V,k);  
    
    if strcmp(sunsky,'sun')
        [W,H,phi,H_shifted] = sun_solve(V,C,W,H);
    elseif strcmp(sunsky,'sky')
        [n,f] = size(V);
        phi = zeros(n,1);
        H_shifted = zeros(f,n);
        [W,H] = sky_solve(V,C,W,H);
    else
        ME = MException('Invalid input for sunsky',sunsky);
        throw(ME)
    end
end
function [W,H] = sky_solve(V,C,W,H)
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
    figure 
    subplot(2,2,1);
    frames =1:m;
    plot(frames,V,'r');
    hold
    plot(frames,Hbasis,'b');
    plot(frames,H,'c');
    xlabel('time');
    ylabel('intensity');
    title('comparing our ACLS to toolbox');
    legend('pixel intensity','their H curve','our H curve');
    %clear NewA    
end

function [W_f,H_f,phi, H_shifted] = sun_solve(V,C,W,H)
    %n = # pixels ;  m = # frames 
    [n,m] = size(V);
    %initialize phi, our shift map if sunsky = 'sun'    
     phi = calculatePhi(V,H,n,m);
    
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
            V_shift = shift_V(V, m, phi); 
            %{
            V_shift = zeros(n,1);
            for d = 1:n 
                %new column to grab data from(clamped between 1 and m) 
                k = max(min((i + phi(d)),m),1);
                V_shift(d) = V(d,k);
            end
            %}
           
            %to include confidence, M = C_i * W and d = C_i * V_i
            M = double(C(:,i) .* W);
            d = double(C(:,i) .* V_shift);
            
            x = lsqnonneg(M,d);
            H(:,i) = x;
        end
        
        fprintf('optimizing phi \n');
        %now we solve for phi
         phi = calculatePhi(V,H,n,m);
         
         fprintf('creating shifted H matrix \n');
         %return full shifted H matrix
         H_shifted = getShiftedMatrix(H,phi);
        
    end
    
    i = round(n / 2);
    figure
    frames = 1:m;
    colors = ['c','y','g','b'];
    %shifted_frames = frames + phi(i,:);
        hold on;
        
    for it = 1:2
        shifted_frames = frames - phi(i,:);

        phi(i,:);
        plot(frames, C(i,:) .* V(i,:),colors(2 + it));
        
        plot(shifted_frames, C(i,:) .* V(i,:),colors(it));
        i = i + round(n/8);
    end
            plot(frames, H,'r');
            
            
     i = round(n/3) + 7;
     figure
     for it = 1:4
            frames = 1:m; 
            shift_frames = (1 + phi(i)):(m + phi(i));
            
            subplot(2,2,it)
            hold on;
            plot(frames, (C(i,:) .* V(i,:)),'m');
            plot(frames, (W(i,:) * H),'r');
            plot(shift_frames, (W(i,:) * H),'c');
            title(strcat('Pixel #',num2str(i),'tracked'));
            xlabel('time');
            ylabel('pixel value');
            legend('pixel intensities','sun illumination','shifted curve');
            i = i + round(n/8);
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

function phi_n = calculatePhi(V,H,n,f)
phi_n = zeros(n,1);
    for i = 1:n
        [r,lag] = xcorr(V(i,:), H);
        [~,I] = max(abs(r));
        lagDiff = lag(I);
        phi_n(i) = lagDiff;
    end
    %{
    frames = 1:f;
    figure
    plot(frames, V(10,:), 'r');
    hold on;
    plot(frames, H, 'g');
    plot(frames, shiftByPhi(H,phi_n(10),f), 'b');
    phi_n(10)
    %}

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

function H_shifted = getShiftedMatrix(H,phis)
    f = numel(H);
    n = numel(phis);
    H_shifted = zeros(f,n);
    for i = 1:n
        H_shifted(:,i) = shiftByPhi(H,phis(i,:),f);
    end
end

function H_shift = shiftByPhi(H,phi,f)
H_shift = zeros(1,f);
shifted_indices = max(1 + phi , 1) : min(f + phi, f);
H_shift(shifted_indices) = H(:,1:numel(shifted_indices));
end

function V_shift = shift_V(V,j,phi)
    [n,f] = size(V);
    indices = 0:(n-1);
    indices = ((indices * f) + j) + phi';
    for i = 1:n
        indices(i) = clamp(indices(i),(f*i-f+1),(f * i)); 
    end
    V_shift = V(indices)';
end

function clamped = clamp(A,l,u)
    clamped = min(max(A,l),u);
end