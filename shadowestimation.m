function [S, threshold] = shadowestimation(images, skymask)

%   output: f x n matrix of binary values.
%       0 -> pixel is in shadow at this frame.
%       1 -> pixel is in sunlight at this frame
%   threshs: 1 x n matrix of thresholds for each pixel

    k = 1.5; %value from the paper
    p = 0.25; %percent of pixels to count as "shadowed"

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
    threshold = m_min .* k;
    
    
    for i = 1:f
        S(i,:) = images(i,:) > threshold;
    end
    
    % perform bilateral filtering
    %create a m x n x f matrix S_mov
    S_mov = replace_sky(S,skymask);
    
    for i = 1:f
        im = S_mov(:,:,i);

        im = bwmorph(im, 'spur');
        im = bwmorph(im, 'clean');
        im = bwmorph(im, 'fill');

        S_mov(:,:,i) = im;
    end
    
    %display stuff
    frames = 1:f; 
    index = 15000;
    figure

    for it=1:4
        pixel = images(:,index);

        shadow = 50 * S(:,index);
        thresh = threshold(:,index);

        subplot(2,2,it)
        plot(frames, shadow, 'm')
        hold on;
        plot(frames, pixel,'r')
        plot(frames, thresh,'c')
        
        title(strcat('Pixel #',index,' tracked'));
        xlabel('time');
        ylabel('pixel value');
        hold off;
        index=index+25000;
    end

end

function withSky = replace_sky(new_values,sky_mask)
[n,m] = size(sky_mask);
f = size(new_values, 1);
withSky = zeros(n,m,f);
land_indices = find(sky_mask);
    for i = 1:f
    sky_mask(land_indices) = new_values(i,:);
    withSky(:,:,i) = sky_mask;
    end
end
