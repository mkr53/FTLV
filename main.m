function [ output_args ] = main( input_args )

%load matrix of image information
if ~exist('imagematrix.mat')
    fprintf('No image matrix data found. Creating Matrix...')
    [images, image_size] = matrixGenerate;
    save('imagematrix.mat', 'images','image_size');
else
    load('imagematrix.mat', 'images','image_size');
end

%this is our binary sky mask. we will remove sky pixels from the matrix
%that we do calcualtions on 
sky_mask = imresize(imread('skymask.png', 'png'),image_size);

no_sky_images = removeSky(images,sky_mask);

images_orig = images;
images = no_sky_images;

%f = # of frames
%n = # of pixels
[f,n] = size(images);

%compute binary shadow estimation for each frame
size(images(:,:,1))
greyscale_images = rgb2gray(images);
S = shadowestimation(greyscale_images, sky_mask);


A = double(images(:,:,1)) .* (1-S);
%factorize F(t) into Sky
[W_sky, H_sky] = nnmf(A, 1);

%to display the image, we must replace the sky
frame = replaceSky(images, S, sky_mask);

index = 110000;
displayAppearance(images, S, index);


implay(frame)

end

function non_sky =  removeSky(images, sky_mask)
%return new matrix: images matrix with all sky pixels removed
%change sky_mask so that it is a single row
mask = reshape(sky_mask, 1, numel(sky_mask));
land_indices = find(mask);
non_sky = images(:,land_indices,:);
end

function withSky = replaceSky(images, new_values, sky_mask)
[n,m] = size(sky_mask);
f = size(new_values, 1);
withSky = zeros(n,m,f);
land_indices = find(sky_mask);
    for i = 1:f
    sky_mask(land_indices) = new_values(i,:);
    withSky(:,:,i) = sky_mask;
    end
end