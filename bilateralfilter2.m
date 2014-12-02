function [ filteredimage] = bilateralfilter2( image ,winsize)
%apply a bilateral filter to a binary image- image and return the filtered
%image
[a,b]=size(image);
filteredimage=image;
for i = winsize+1:(a-winsize),
    for j= winsize+1:(b-winsize),
        g=0;
        for w=1:winsize;
        deln=image(i+w,j)-image(i,j);
        deln=deln*deln;
        dels=image(i-w,j)-image(i,j);
        dels=dels*dels;
        delw=image(i,j-w)-image(i,j);
        delw=delw*delw;
        dele=image(i,j+w)-image(i,j);
        dele=dele*dele;
        g=g+image(i,j)+deln+dels+delw+dele;
        end
        g=round(g/winsize*5);
        filteredimage(i,j)=g;
    end
end
end

