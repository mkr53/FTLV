%VideoReader.getFileFormats()

workingDir = '/Users/elisecotton/Documents/MATLAB/FTLV/data/arizona/';

if ~exist(strcat(workingDir,'images'), 'dir')
    mkdir(workingDir)
    mkdir(workingDir,'images')
end

%shuttleVideo = VideoReader('shuttle.avi');
shuttleVideo = VideoReader('arizonadata.mpg');

numberOfFrames = shuttleVideo.NumberOfFrames;
shuttleVideo.duration
%hasFrame(shuttleVideo)

%grab height and width of images
num_pixels = shuttleVideo.height * shuttleVideo.width;

%our main images matrix


ii = 1;

while ii < numberOfFrames
    img = read(shuttleVideo,ii);

    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
    imwrite(img,fullfile(workingDir,'images',sprintf('img%d.jpg',ii)));
    
    ii = ii+1;
end

imageNames = dir(fullfile(workingDir,'images','*.jpg'));
imageNames = {imageNames.name}';

outputVideo = VideoWriter(fullfile(workingDir,'shuttle_out.avi'));
outputVideo.FrameRate = shuttleVideo.FrameRate;
open(outputVideo)

for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,'images',imageNames{ii}));
   writeVideo(outputVideo,img)
end

close(outputVideo)
