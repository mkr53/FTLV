function [ output_args ] = main( input_args )
%MAIN Summary of this function goes here
%   Detailed explanation goes here
a=input_args;
%set up temp working folder for image sequence
workingDir = tempname;
%create videoreader
mkdir(workingDir)
mkdir(workingDir,'images')
shuttleVideo = VideoReader(a);
%create image sequence
ii = 1;
while hasFrame(shuttleVideo)
   img = readFrame(shuttleVideo);
   filename = [sprintf('%03d',ii) '.jpg'];
   fullname = fullfile(workingDir,'images',filename);
   imwrite(img,fullname)    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
   ii = ii+1;
end
%set up file names in matrice so we can access later
imageNames = dir(fullfile(workingDir,'images','*.jpg'));
%imageNames is a Nx1 cell matrice of file names
imageNames = {imageNames.name}';
luminanceperpixel=shadowestimation(imageNames);
end

