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
    clear images

    %run FTLV algorithm
    [I_sky_r, W_sky_r, H_sky_r, I_sun_r, S_r, W_sun_r, H_sun_r, phi_r, threshs_r, H_shifted_r] = FTLV(F_t(:,:,1),'r');
    [I_sky_g, W_sky_g, H_sky_g, I_sun_g, S_g, W_sun_g, H_sun_g, phi_g, threshs_g, H_shifted_g] = FTLV(F_t(:,:,2),'g');
    [I_sky_b, W_sky_b, H_sky_b, I_sun_b, S_b, W_sun_b, H_sun_b, phi_b, threshs_b, H_shifted_b] = FTLV(F_t(:,:,3),'b');
    
   
    %total reconstracted matrix = 
    F = I_sky + S_r*I_sun;
    %RMS error
    RMS= calculateRMS(F, F_t);
    disp(strcat('RMS: ',num2str(RMS)));
    
    %normalize phi for visualizing
    if min(phi_r) < 0
        phi_r = phi_r - min(phi_r(:));
    end
    phi_r = (phi_r - min(phi_r(:))) ./ (max(phi_r(:)) - min(phi_r(:)));
    phi_r = imcomplement(phi_r);

    %plot sun illumination
    frame_r = W_sun_r;
    frame_r = backIntoImage(frame_r, sky_mask);
    figure
    subplot(3,1,1)
    imshow(frame_r);
    title('W_sun_r=sun illumination image');
    
    phis_r = phi_r;
    phis_r = backIntoImage(phis_r, sky_mask);
    subplot(3,1,2);
    imshow(phis_r);
    title('shift map');
    
    frames = 1:f;
    subplot(3,1,3);
    plot(frames, H_sun_r, 'c');
    xlabel('time');
    ylabel('sun illumination');
    title('sun illumination over time');
    
    %plot sky illumination
    frame_r = W_sky_r;
    frame_r = backIntoImage(frame_r, sky_mask);
    figure
    subplot(3,1,1)
    imshow(frame_r);
    title('W_sky_r=sky illumination image');
    
    frames = 1:f;
    subplot(3,1,2);
    plot(frames, H_sky_r, 'c');
    xlabel('time');
    ylabel('sky illumination');
    title('sky illumination over time');
    %{
    S_mov = replaceSky(S_r, sky_mask);
    figure
    implay(S_mov)
    
    sky_frame = backIntoImage(W_sky_r,sky_mask);
    figure
    imshow(sky_frame)
    %}
    figure
    imshow(backIntoImage(W_sun_r, sky_mask));
    figure
    imshow(backIntoImage(W_sun_g, sky_mask));
    figure
    imshow(backIntoImage(W_sun_b, sky_mask));
    
    color_sky = createColorImage(W_sky_r, W_sky_g, W_sky_b, sky_mask);
    color_sun = createColorImage(W_sun_r, W_sun_g, W_sun_b, sky_mask);
    
    figure 
    imshow(color_sky)
    figure
    imshow(color_sun)
    %displayAppearance(F_t, S_sun_r, threshs, I_sky',I_sun, 'intensity');

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

function withSky = replaceSky(new_values, sky_mask)
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

function color_image = createColorImage(R,G,B,sky_mask)
%to display the image, we must replace the sky
    image_size = size(sky_mask);
    frame_r = (R - min(R(:))) ./ (max(R(:)) - min(R(:)) );
    frame_r = backIntoImage(frame_r, sky_mask);
    frame_g = G;
    frame_g = (frame_g - min(frame_g(:))) ./ (max(frame_g(:)) - min(frame_g(:)) );
    frame_g = backIntoImage(frame_g, sky_mask);
    frame_b = B;
    frame_b = (frame_b - min(frame_b(:))) ./ (max(frame_b(:)) - min(frame_b(:)) );
    frame_b = backIntoImage(frame_b, sky_mask);
    color_image = zeros(image_size(1),image_size(2), 3);
    color_image(:,:,1) = frame_r;
    color_image(:,:,2) = frame_g;
    color_image(:,:,3) = frame_b;
    color_image = color_image * 256; 
end