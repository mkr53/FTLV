function [ totMSE] = calculateRMS( theirinp,ourout )
%return RMS error given original input and calculated after decomposition
[frames,pixels,color]=size(theirinp);
size(theirinp)
MSE=zeros(pixels,1);
for pix=1:pixels,
    ourframe=ourout(:,pix,1);
    theirframe=theirinp(:,pix,1);
    D = abs(ourframe-theirframe).^2;
    MSE(pix,1) = sqrt(sum(D(:))/frames);
end
totMSE = sum(MSE(:,1))/pixels;
end
