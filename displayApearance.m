data_directory = '/Users/elisecotton/Documents/MATLAB/FTLV/data';

if ~exist([data_directory 'imagematrix.mat'])
    images = matrixGenerate;
    save('imagematrix.mat', 'images');
else
    images = load('imagematrix.mat', 'images');
end

num_frames = size(images, 1);

frames = 1:num_frames; 

pixel_200_red = images(:,200,1);

figure
plot(pixel_200_red, frames)