function displayAppearance(images, shadows , threshs, sky, sun, ylabeltitle)

%images

num_frames = size(images, 1);

frames = 1:num_frames; 
%index=110000;
index = 15000;
figure

for it=1:4,
pixel = images(:,index);


shadow_200 = 50 * shadows(:,index);
thresh = threshs(:,index);
sky_index = sky(:,index);

subplot(2,2,it)
plot(frames, shadow_200, 'm')

hold on;
plot(frames, pixel,'r')
plot(frames, thresh,'c')
title(strcat('Pixel #',index,' tracked'));
xlabel('time');
ylabel(ylabeltitle);
hold off;
index=index+25000;
end

index = 15000;
figure
for it=1:4,
pixel = images(:,index);

sky_index = sky(:,index);
sun_index = sun(:,index);
thresh = threshs(:,index);


subplot(2,2,it)
plot(frames, pixel, 'r')

hold on;
plot(frames, thresh, 'b');
plot(frames, sky,'c')
plot(frames, sun,'y')
title(strcat('Pixel #',index,' tracked'));
xlabel('time');
ylabel(ylabeltitle);
hold off;
index=index+25000;
end



end