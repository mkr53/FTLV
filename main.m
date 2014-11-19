function [ output_args ] = main( input_args )

%load matrix of image information
if ~exist('imagematrix.mat')
    fprintf('No image matrix data found. Creating Matrix...')
    [images, image_size] = matrixGenerate;
    save('imagematrix.mat', 'images','image_size');
else
    load('imagematrix.mat', 'images','image_size');
end

[f,n] = size(images);
%compute binary shadow estimation for each frame
size(images(:,:,1))
S = shadowestimation(images(:,:,1));

halfway = floor(0.5 * f);
frame = reshape(S(halfway,:), image_size);

figure
imshow(frame)

end

