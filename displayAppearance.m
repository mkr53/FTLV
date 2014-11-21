function displayAppearance(images, shadows , sky, ylabeltitle)

%images

num_frames = size(images, 1);

frames = 1:num_frames; 
%index=110000;
index = 100;
figure

for it=1:4,
pixel_200_red = rgb2gray(images(:,index,:));
shadow_200 = 50 * shadows(:,index);
sky_index = sky(:,index);

subplot(2,2,it)
plot(frames, pixel_200_red,'r')
hold on;
plot(frames, shadow_200, 'm')
title('Pixel #%d tracked');
xlabel('time');
ylabel(ylabeltitle);
plot(frames, sky_index, 'b')
hold off;
index=index+25000;
end



end