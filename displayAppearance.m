function displayAppearance(images, shadows, index)

%images

num_frames = size(images, 1);

frames = 1:num_frames; 


pixel_200_red = rgb2gray(images(:,index,:));
shadow_200 = 50 * shadows(:,index);

figure
plot(frames, pixel_200_red,'r')
hold
plot(frames, shadow_200, 'm')

end