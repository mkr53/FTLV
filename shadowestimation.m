function [ luminanceperpixel ] = shadowestimation( imagefiles )

%load matrix of image information
if ~exist('imagematrix.mat')
    images = matrixGenerate;
    save('imagematrix.mat', 'images');
else
    load('imagematrix.mat', 'images');
end

%UNTITLED Summary of this function goes here
%   output: f x n matrix of binary values. 
%       0 -> pixel is in shadow at this frame. 
%       1 -> pixel is in sunlight at this frame
b=imread(char(imagesfiles(1)));
[numr,numc]=size(b);
luminanceperpixel=zeros(length(imagefiles),numr*numc);
for i=1:length(imagefiles),
    %cim is rowxcolumnx3
    cim=imread(char(imagefiles(i)));
    rim=cim(:,:,1);
    vim=imreshape(rim,1,:);
    luminanceperpixel(i,:)=vim;
        

end
j=[1:length(imagefiles)];
plot(j,luminanceperpixel(j,1));

end
