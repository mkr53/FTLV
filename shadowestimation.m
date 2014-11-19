function [S] = shadowestimation()

%UNTITLED Summary of this function goes here
%   output: f x n matrix of binary values. 
%       0 -> pixel is in shadow at this frame. 
%       1 -> pixel is in sunlight at this frame
k = 1.5; %value from the paper
p = 0.2; %percent of pixels to count as "shadowed"

images = randi([0 256], 5, 25);
images

[f,n] = size(images);
num_shadowed = floor(f*p)

%sort each column (pixel data) 
[sorted,I] = sort(images,1);
sorted

%CURRENT ISSUE: indexing in I is only relative to the column that its in..
%how to make it correct in terms of larger matrix
shadowed_pixels = I(1:num_shadowed,:)

S = zeros(f,n);
S(shadowed_pixels) = 1;

end
