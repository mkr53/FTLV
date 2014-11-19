data_directory = '/Users/elisecotton/Documents/MATLAB/FTLV/data';

exist('imagematrix.mat')

if ~exist('imagematrix.mat')
    [images, image_size]= matrixGenerate;
    save('imagematrix.mat', 'images', 'image_size');
else
    load('imagematrix.mat', 'images');
end
%images
size(images)

num_frames = size(images, 1);

frames = 1:num_frames; 

pixel_200_red = images(:,200,1);

figure
plot(frames, pixel_200_red)