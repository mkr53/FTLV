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
    
    
    for i = 1:f
        S(i,:) = images(i,:) > threshold;
    end

end
