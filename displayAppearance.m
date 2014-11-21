function displayAppearance(images, shadows, sky, index)

%images

num_frames = size(images, 1);

frames = 1:num_frames; 


pixel_200_red = (images(:,index,1));
shadow_200 = 250 * shadows(:,index);
%shadow_200 = shadows(:,index);

sky_index = sky(:,index); 


figure
plot(frames, pixel_200_red,'r')
hold on
plot(frames, shadow_200, 'm')

plot(frames, sky_index, 'b')
hold off

end