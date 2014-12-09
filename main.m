function main
    %load matrix of image information
    %images will be a f x (m x n) x 3 matrix where f is the number of frames and
    %mxn is the size of an image
    
    %NOTE: to load new data set, you must delete imagematrix.mat
    if ~exist('imagematrix.mat')
        fprintf('No image matrix data found. Creating Matrix... \n')
        [images, image_size] = matrixGenerate;
        save('imagematrix.mat', 'images','image_size');
    else
        fprintf('loading image matrix \n');
        load('imagematrix.mat', 'images','image_size');
    end

    %this is our binary sky mask. we will remove sky pixels from the matrix
    %that we do calcualtions on 
    sky_mask = imresize(imread('skymask.png', 'png'),image_size);

    %F_t is a fxnx3 matrix where f is the number of frames and n is the number
    %of NON-SKY pixels
    F_t = removeSky(images,sky_mask);
    [f,n,d] = size(F_t);

    %run FTLV algorithm
    [I_sky, W_sky_r, H_sky_r, I_sun, S_sun_r, W_sun_r, H_sun_r, phi_r, threshs] = FTLV(F_t(:,:,1), sky_mask);
    %[W_sky_g, H_sky_g, S_sun_g, W_sun_g, H_sun_g, phi_g] = FTLV(F_t(:,:,2), sky_mask);
    %[W_sky_b, H_sky_b, S_sun_b, W_sun_b, H_sun_b, phi_b] = FTLV(F_t(:,:,3), sky_mask);
    
    %normalize phi for visualizing
    if min(phi_r) < 0
        phi_r = phi_r - min(phi_r(:));
    end
    phi_r = (phi_r - min(phi_r(:))) ./ (max(phi_r(:)) - min(phi_r(:)));

    
    frame_r = W_sun_r;
    frame_r = backIntoImage(frame_r, sky_mask);
    figure
    subplot(3,1,1)
    imshow(frame_r);
    
    phis_r = phi_r;
    phis_r = backIntoImage(phis_r, sky_mask);
    subplot(3,1,2);
    imshow(phis_r);
    
    frames = 1:f;
    subplot(3,1,3);
    plot(frames, H_sun_r, 'c');
    %{
    %to display the image, we must replace the sky
    frame_r = W_sky_r;
    frame_r = (frame_r - min(frame_r(:))) ./ (max(frame_r(:)) - min(frame_r(:)) );
    frame_r = backIntoImage(frame_r, sky_mask);
    frame_g = W_sky_g;
    frame_g = (frame_g - min(frame_g(:))) ./ (max(frame_g(:)) - min(frame_g(:)) );
    frame_g = backIntoImage(frame_g, sky_mask);
    frame_b = W_sky_b;
    frame_b = (frame_b - min(frame_b(:))) ./ (max(frame_b(:)) - min(frame_b(:)) );
    frame_b = backIntoImage(frame_b, sky_mask);
    frame = zeros(image_size(1),image_size(2), 3);
    frame(:,:,1) = frame_r;
    frame(:,:,2) = frame_g;
    frame(:,:,3) = frame_b;


    figure
    imshow(frame)
    %}
    displayAppearance(F_t, S_sun_r, threshs, I_sky',I_sun', 'intensity');

end



function non_sky =  removeSky(images, sky_mask)
%return new matrix: images matrix with all sky pixels removed
%change sky_mask so that it is a single row
mask = reshape(sky_mask, 1, numel(sky_mask));
land_indices = find(mask);
non_sky = images(:,land_indices,:);
end

function new_S = mov_to_matrix(S_mov, sky_mask)
S_mov = permute(S_mov,[3 1 2]);
[f,n,m] = size(S_mov);
S_mov = reshape(S_mov, f, (n*m));
mask = reshape(sky_mask, 1, numel(sky_mask));

land_indices = find(mask);
new_S = round(S_mov(:,land_indices));

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