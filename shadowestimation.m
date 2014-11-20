function [S] = shadowestimation(images)

%UNTITLED Summary of this function goes here
%   output: f x n matrix of binary values. 
%       0 -> pixel is in shadow at this frame. 
%       1 -> pixel is in sunlight at this frame
k = 1.5; %value from the paper
p = 0.2; %percent of pixels to count as "shadowed"

%images = randi([0 256], 5, 25);

[f,n] = size(images);
num_shadowed = floor(f*p);

%sort each column (pixel data) 
[sorted,I] = sort(images,1);

%first, compute the median value m_min of the num_shadowed smallest intensities at
%each pixel
shadowed_pixels = sorted(1:num_shadowed,:);
m_min = median(shadowed_pixels, 1);

%set shadow function to one for each F(t) > k*m_min and zero otherwise
S = zeros(f,n);

threshold = m_min * k; 
size(threshold)
size(images)
for i = 1:f
    S(i,:) = images(i,:) > threshold; 
end


%{
columns = 1:(n*num_shadowed);
columns = ceil(columns / num_shadowed);

shadowed_pixels = reshape(shadowed_pixels, 1, numel(shadowed_pixels));
linearInd = sub2ind([f n], shadowed_pixels, columns);
S = ones(f,n);
S(linearInd) = 0;
%}

end
