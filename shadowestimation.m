function [S] = shadowestimation()

%load matrix of image information
if ~exist('imagematrix.mat')
    images = matrixGenerate;
    save('imagematrix.mat', 'images');
else
    load('imagematrix.mat', 'images');
end

%UNTITLED Summary of this function goes here
%   output: f x n matrix of binary values. 
%       0 -> pixel is in shadow at this frame. 
%       1 -> pixel is in sunlight at this frame
k = 1.5 %value from the paper
p = 0.2 %percent of pixels to count as "shadowed"

[f,n] = size(images);
%sort each column (pixel data) 
[sorted,I] = sort(images,1);
shadowed_pixels = I(1:(f*p),:);

S = zeros(f,n);
S(shadowed_pixels) = 1;

end
