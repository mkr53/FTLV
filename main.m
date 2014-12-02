function main

%load matrix of image information
if ~exist('imagematrix.mat')
    fprintf('No image matrix data found. Creating Matrix...')
    [images, image_size] = matrixGenerate;
    save('imagematrix.mat', 'images','image_size');
else
    load('imagematrix.mat', 'images','image_size');
end

%images will be a f x (m x n) x 3 matrix where f is the number of frames and
%mxn is the size of an image

%this is our binary sky mask. we will remove sky pixels from the matrix
%that we do calcualtions on 
sky_mask = imresize(imread('skymask.png', 'png'),image_size);

%F_t is a fxnx3 matrix where f is the number of frames and n is the number
%of NON-SKY pixels
F_t = removeSky(images,sky_mask);

%f = # of frames
%n = # of pixels (non-sky)
[f,n] = size(F_t);

%compute binary shadow estimation for each frame
%we will do this individualy for the red, green, and blue components
fprintf('estimating shadows \n');
greyscale_images = rgb2gray(F_t);

[S,threshs] = shadowestimation(greyscale_images);
[S_r, threshs_r] = shadowestimation(F_t(:,:,1));
[S_g, threshs_g] = shadowestimation(F_t(:,:,2));
[S_b, threshs_b] = shadowestimation(F_t(:,:,3));


S_mov = replaceSky(F_t, S, sky_mask);

fprintf('applying bilateral filter \n');
%here is where we will perform the bilateral filter
%{
for i = 1:f
im = S_mov(:,:,f);
im = bilateralfilter2(im, 1);
S_mov(:,:,f) = im;
end
%}
implay(S_mov)

%S_mov = reshape(S_mov, image_size(1)*image_size(2), f);
%S_mov = permute(S_mov, [2 1]);

%S = removeSky(S_mov, sky_mask);


%factorize F(t) into Sky: W_sky and H_sky
fprintf('finding I_sky \n');
A = double(F_t(:,:,1)');
[W_sky_r,H_sky_r] = ACLS(A, (1-S)', 'sky');

I_sky = W_sky_r * H_sky_r;


%next, we factorize F(t) - I(sky) = I(sun) into its W_sun and H_sun parts
%(clamped to 0)

I_sun = max(double(F_t(:,:,1)) - double(I_sky'), 0);
fprintf('finding I_sun \n');
A = double(F_t(:,:,1)');
[W_sun_r, H_sun_r] = ACLS(A, S', 'sun');

%to display the image, we must replace the sky
%frame = replaceSky(images, 256 * H_sky, sky_mask);
frame_r = W_sky_r;
frame_r = (frame_r - min(frame_r(:))) ./ (max(frame_r(:)) - min(frame_r(:)) );
frame = backIntoImage(frame_r, sky_mask);



%displayAppearance(F_t, S, I_sky', index);
displayAppearance(F_t, S_r, threshs_r, I_sky', 'intensity');

%frame
figure
imshow(frame)
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