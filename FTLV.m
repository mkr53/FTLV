function [I_sky, W_sky, H_sky, I_sun, S, W_sun, H_sun, phi, threshs] = FTLV(F_t, sky_mask)

%f = # of frames
%n = # of pixels (non-sky)
[f,n] = size(F_t);

%initialize outputs
W_sky = zeros(n,1);
H_sky = zeros(1,f);
S_sun = zeros(n,f);
W_sun = zeros(n,1);
H_sun = zeros(1,f);
phi = zeros(n,1);

%compute binary shadow mask S
fprintf('estimating shadows \n');
[S, threshs] = shadowestimation(F_t, sky_mask);


%S_mov = replaceSky(F_t, S, sky_mask);
%{
fprintf('applying bilateral filter \n');
if ~exist('filteredshadow.mat')
    fprintf('No filtered shadow data found. Creating Matrix... \n')
    fprintf('applying bilateral filter \n');
    w     = 5;       % bilateral filter half-width
    sigma = [3 0.1]; % bilateral filter standard deviations
    for i = 1:f
        
        im = S_mov(:,:,i);
      
        im = bfilter2(im, w, sigma);
        gim=round(im);
        S_mov(:,:,i) = gim;
    end

    save('filteredshadow.mat', 'S_mov');
else
    fprintf('loading filtered shadow \n');
    load('filteredshadow.mat', 'S_mov');
end
%}

%here is where we will perform the bilateral filter
% Set bilateral filter parameters.
%w     = 5;       % bilateral filter half-width
%sigma = [3 0.1]; % bilateral filter standard deviations
%{
for i = 1:f
    im = S_mov(:,:,i);

    im = bwmorph(im, 'spur');
    im = bwmorph(im, 'clean');
    im = bwmorph(im, 'fill');

    S_mov(:,:,i) = im;
end

%}
%S = mov_to_matrix(S_mov, sky_mask);


%S = removeSky(S_mov, sky_mask);
    %A = double(F_t');
%{
if ~exist('skyest.mat')
    fprintf('No sky data found. Creating Matrix... \n')
    %factorize F(t) into Sky: W_sky and H_sky
    fprintf('finding I_sky \n');
    A = double(F_t');
    [W_sky,H_sky] = ACLS(A, (1-S)', 'sky');

    I_sky = W_sky * H_sky;
    save('skyest.mat', 'W_sky','H_sky','I_sky');
else
    fprintf('Loading sky data \n');
    load('skyest.mat', 'W_sky','H_sky','I_sky');
end
%}
    %factorize F(t) into Sky: W_sky and H_sky
    fprintf('finding I_sky \n');
    A = double(F_t');
    [W_sky,H_sky] = ACLS(A, (1-S)', 'sky');

    I_sky = W_sky * H_sky;
%{
%use this version, to compare decomposition 
NewA=A .* (1-S)';
fprintf('size of input matrix: \n');
size(NewA)
[Wcoeff,Hbasis,numIter,tElapsed,finalResidual]=wnmfrule(NewA,1);
%}


%next, we factorize F(t) - I(sky) = I(sun) into its W_sun and H_sun parts
I_sun = max(double(F_t) - double(I_sky'), 0);
fprintf('finding I_sun \n');
%A = double(F_t');
%[W_sun_r, H_sun_r, phi] = ACLS(I_sun', S', 'sun');



end