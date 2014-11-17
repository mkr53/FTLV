function [ luminanceperpixel ] = shadowestimation( imagefiles )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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
