function displayAppearance(images, shadows , threshs, sky, ylabeltitle)

%images

num_frames = size(images, 1);

frames = 1:num_frames; 
%index=110000;
index = 15000;
figure

for it=1:4,
pixel_200_red = images(:,index,1);
pixel_200_green = images(:,index,2);
pixel_200_blue = images(:,index,3);

shadow_200 = 50 * shadows(:,index);
thresh = threshs(:,index);
sky_index = sky(:,index);

subplot(2,2,it)
plot(frames, shadow_200, 'm')

hold on;
plot(frames, pixel_200_red,'r')
%plot(frames, pixel_200_green,'g')
%plot(frames, pixel_200_blue,'b')

plot(frames, thresh,'c')
title(strcat('Pixel #',index,' tracked'));
xlabel('time');
ylabel(ylabeltitle);
%plot(frames, sky_index, 'y')
hold off;
index=index+25000;
end



end