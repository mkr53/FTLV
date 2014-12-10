function [I_sky, W_sky, H_sky, I_sun, S, W_sun, H_sun, phi, threshs] = FTLV(F_t, sky_mask)

%f = # of frames
%n = # of pixels (non-sky)
[f,n] = size(F_t);

%initialize outputs
I_sky = zeros(n,f);
I_sun = zeros(n,f);
W_sky = zeros(n,1);
H_sky = zeros(1,f);
W_sun = zeros(n,1);
H_sun = zeros(1,f);
phi = zeros(n,1);

if ~exist('shadowest.mat')
        fprintf('No shadow data found. Creating Matrix... \n')
        fprintf('estimating shadows \n');
        [S, threshs] = shadowestimation(F_t);
        save('shadowest.mat', 'S','threshs');
    else
        fprintf('loading shadow matrix \n');
        load('shadowest.mat', 'S','threshs');
    end
%compute binary shadow mask S
%fprintf('estimating shadows \n');
%[S, threshs] = shadowestimation(F_t, sky_mask);


%factorize F(t) into Sky: W_sky and H_sky

if ~exist('skyest.mat')
        fprintf('No sky data found. Factorizing Matrix... \n')
        fprintf('finding I_sky \n');
        A = double(F_t');
        [W_sky,H_sky] = ACLS(A, (1-S)', 'sky');

        I_sky = W_sky * H_sky;
        save('skyest.mat', 'W_sky','H_sky','I_sky');
    else
        fprintf('loading sky decomposition \n');
        load('skyest.mat', 'W_sky','H_sky','I_sky');
end
%{    
fprintf('finding I_sky \n');
A = double(F_t');
[W_sky,H_sky] = ACLS(A, (1-S)', 'sky');

I_sky = W_sky * H_sky;
%}
%{
%use this version, to compare decomposition 
NewA=A .* (1-S)';
fprintf('size of input matrix: \n');
size(NewA)
 
[Wcoeff,Hbasis,numIter,tElapsed,finalResidual]=wnmfrule(NewA,1);
%}
  

%next, we factorize F(t) - I(sky) = I(sun) into its W_sun and H_sun parts
I_sun = max(double(F_t) - double(I_sky'), 0);

if ~exist('sunest.mat')
        fprintf('No sun data found. Factorizing Matrix... \n')
        fprintf('finding I_sun \n');
        [W_sun, H_sun, phi] = ACLS(I_sun', S', 'sun');

        save('sunest.mat', 'W_sun','H_sun', 'phi');
    else
        fprintf('loading sun decomposition \n');
        load('sunest.mat', 'W_sun','H_sun', 'phi');
end
    
    fprintf('FTLV complete \n');
   %construct the I_sun matrix (nested for loops because shift map makes working with
        %matrices complicated)
        %THIS TAKES A REALLY LONG TIME... FIX IT?
      %{  
          for i = 1:n
            i
            w = W_sun(i,:);
            for j = 1:f
                j
                shifted_index = min(max(j + phi(i,:),1),f);
                I_sun(j,i) = w * H_sun(:,shifted_index);
            end
      

            end
         %}
%{
I_sun = max(double(F_t) - double(I_sky'), 0);
fprintf('finding I_sun \n');
[W_sun, H_sun, phi] = ACLS(I_sun', S', 'sun');
%}

end

