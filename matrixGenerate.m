% this function gives us images an fxnx3 matrix where f is the number of
% frames and n is the number of pixels in a single frame
function images = matrixGenerate
    
dataVideo = VideoReader('arizonadata.mpg');

num_frames = dataVideo.NumberOfFrames;
num_pixels = dataVideo.height * dataVideo.width;

%our main images matrix
%images = zeros(num_frames, num_pixels,3);

matrix = reshape(read(dataVideo),[num_pixels 3 num_frames]);
matrix = permute(matrix,[3 1 2]);

size(matrix)

images = matrix;

end

