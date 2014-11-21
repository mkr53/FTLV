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

fprintf('estimating shadows')
%compute binary shadow estimation for each frame
greyscale_images = rgb2gray(images);
S = shadowestimation(greyscale_images);
size(S)

fprintf('factorizing F')
A = double(images(:,:,1));
%factorize F(t) into Sky
[W_sky_r,H_sky_r] = ACLS(A, 1-S);

A = double(images(:,:,2));
[W_sky_g,H_sky_g] = ACLS(A, 1-S);

A = double(images(:,:,3));
[W_sky_b,H_sky_b] = ACLS(A, 1-S);


H_sky = zeros(image_size(1),image_size(2), 3);





%to display the image, we must replace the sky
%frame = replaceSky(images, 256 * H_sky, sky_mask);
frame_r = backIntoImage(H_sky_r, sky_mask) ./ 256;
frame_g = backIntoImage(H_sky_g, sky_mask) ./ 256;
frame_b = backIntoImage(H_sky_b, sky_mask) ./ 256;

H_sky(:,:,1) = frame_r;
H_sky(:,:,2) = frame_g;
H_sky(:,:,3) = frame_b;

size(H_sky)

%index = 110000;
%displayAppearance(images, S, index);
%frame
imshow(H_sky)
%implay(frame)

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

function withSky = backIntoImage(new_values, sky_mask)
[n,m] = size(sky_mask);
withSky = zeros(n,m);

land_indices = find(sky_mask);
withSky(land_indices) = new_values(:);
end