function [S, threshold] = shadowestimation(images)

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
    
    %smoothed images
    smimages=zeros(f,n);
    
    for i = 1:n
        pixel = images(:,i);
        smimages(:,i) = smooth(double(pixel),61, 'sgolay',1);
    end
    
    for i = 1:f
        %S(i,:) = images(i,:) > threshold;
        S(i,:) = smimages(i,:) > threshold;
    end
    
    %display stuff
    frames = 1:f; 
    index = 15000;
    figure

    for it=1:4

        pixel = images(:,index);
        %y = filter(b,a,double(pixel));


        shadow = 50 * S(:,index);
        thresh = threshold(:,index);

        subplot(2,2,it)
        plot(frames, shadow, 'b')
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
