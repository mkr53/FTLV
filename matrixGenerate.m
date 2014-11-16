
function matrixGenerate
    data_directory = '/Users/elisecotton/Documents/MATLAB/FTLV/data/pyramid/';
    image_directory = [data_directory '2011.05'];

    
    %% read in images.
    disp('  Building image array...')

    ims = zeros(height, width, 3, n, 'uint8');
    index = 1;
    chosen = zeros(n,1);
    for ix=select'
        if index > n
            continue
        end

        name = image_locations(ix).name;
        new_im = imresize(imread(name), [height, width]);

        ims(:,:,:,index) = new_im;
        chosen(index) = ix;
        index = index + 1;
    end

end

