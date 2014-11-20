function [S] = shadowestimation(images, sky_mask)

%UNTITLED Summary of this function goes here
%   output: f x n matrix of binary values. 
%       0 -> pixel is in shadow at this frame. 
%       1 -> pixel is in sunlight at this frame
k = 1.5; %value from the paper
p = 0.25; %percent of pixels to count as "shadowed"

%images = randi([0 256], 5, 25);

[f,n] = size(images);
num_shadowed = floor(f*p);

%sort each column (pixel data) 
sorted = sort(images,1);

%first, compute the median value m_min of the num_shadowed smallest intensities at
%each pixel
shadowed_pixels = sorted(1:num_shadowed,:);
m_min = median(shadowed_pixels, 1);

%set shadow function to one for each F(t) > k*m_min and zero otherwise
S = zeros(f,n);
threshold = m_min * k; 
for i = 1:f
    S(i,:) = images(i,:) > threshold; 
end
for j = 1:n
    h = fspecial('prewitt');    
    B = imfilter(S(:,j)', h);
    S(:,j) = B';
end

%edge-preserving bi-lateral filter
%"Bilateral Filtering for Gray and Color Images" [Tomasi and Manduchi 1998]
%S_f = shadowFilter(S, sky_mask);

%S = S_f;
end


%IGNORE THE NEXT COUPLE FUNCTIONS. THIS FILTER IS HARDER THAN ANTICIPATED
function S_f = shadowFilter(S, sky_mask)
w     = 5;       % bilateral filter half-width
sigma_d = 3;
sigma_r = 0.1;
% Pre-compute Gaussian distance weights.
[X,Y] = meshgrid(-w:w,-w:w);
G = exp(-(X.^2+Y.^2)/(2*sigma_d^2));

%1) reshape S matrix into nxmxf matrix withSky
    [n,m] = size(sky_mask);
    f = size(S, 1);
    withSky = zeros(n,m,f);
    land_indices = find(sky_mask);
    for i = 1:f
        i
        sky_mask(land_indices) = S(i,:);
        %filter image
        h = fspecial('prewitt')
        
        B = imfilter(sky_mask, h);
        %B = bfltGray(double(sky_mask),w,sigma_r,G);
        withSky(:,:,i) = B;
    end
    
    S_f = withSky; 
end

% Implements bilateral filtering for grayscale images.
function B = bfltGray(A,w,sigma_r,G)


% Create waitbar.
%h = waitbar(0,'Applying bilateral filter...');
%set(h,'Name','Bilateral Filter Progress');

% Apply bilateral filter.
dim = size(A);
B = zeros(dim);
for i = 1:dim(1)
   for j = 1:dim(2)
      
         % Extract local region.
         iMin = max(i-w,1);
         iMax = min(i+w,dim(1));
         jMin = max(j-w,1);
         jMax = min(j+w,dim(2));
         I = A(iMin:iMax,jMin:jMax);
      
         % Compute Gaussian intensity weights.
         H = exp(-(I-A(i,j)).^2/(2*sigma_r^2));
      
         % Calculate bilateral filter response.
         F = H.*G((iMin:iMax)-i+w+1,(jMin:jMax)-j+w+1);
         B(i,j) = sum(F(:).*I(:))/sum(F(:));
               
   end
   %waitbar(i/dim(1));
end
end
