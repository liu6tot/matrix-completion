% Include the L1 norm minimizer.
addpath(genpath(pwd));

% Get the image.
image_size = 200;
P = imread('../data/mickey.jpg');
P = imresize(P, [image_size image_size]);
P = im2double(rgb2gray(P));

% Matrix size.
n1 = image_size; n2 = image_size;

% Rank
r = 30;

% Number of elements to observe.
m = round(.60*n1*n2); 
p  = m/(n1*n2);

% Folder the save the results in.
filename = strcat('../results/', num2str(100*p), '_percent_elements/', num2str(r), '_rank/');
mkdir(strcat(filename, '/all_variables/'));
% The file which contains all the errors.
fileID = fopen(strcat(filename, '/result.txt'),'w');

% Generate the actual matrix.
% Generate a low rank approximation of image.
[U,Sigma,V] = svd(full(P), 'econ');
actual_matrix = U(:,1:r)*Sigma(1:r,1:r)*V(:,1:r)';
% Save the image.
% Write the original image.
imwrite(actual_matrix, strcat(filename, '/original_image.png'));

% Test DCT Matrix.
A = DCT_Matrix();
At = A';
dct_z = dct(actual_matrix(:));
estimated = A*dct_z;

% Plot the DCT coefficients of the image to see if its sparse.
dct_actual_matrix = dct2(actual_matrix);
sorted_dct = sort(abs(dct_actual_matrix(:)), 'descend');
figure; plot(sorted_dct);

% Data observed.
Omega_idx = randsample(n1*n2,m);
Omega = zeros(n1, n2);
Omega(Omega_idx) = 1;

fprintf('Matrix completion: %d x %d matrix \n', n1,n2);
fprintf('Rank %d, %.1f%% observation \n', r,100*p);


% Set the initial conditions.
rank_estimated = r;
U_i = eye(rank_estimated);
Y_i = eye(n2);
Y_i = Y_i(1:rank_estimated, :);

% Observed image.
observed_image = Omega.*actual_matrix;

% Write the observed image.
imwrite(observed_image, strcat(filename, '/observed_image.png'));
fprintf('Relative Mean Square Error of observed matrix: %.2f%% \n',...
    (norm(observed_image(:) - actual_matrix(:))/norm(actual_matrix(:)))*100);

% First we'll include the L1-norm
L1_norm = 1;
[X, U, Y, ~, ~] = admm(U_i, Y_i, observed_image, observed_image, Omega, L1_norm);

Z = X*U*Y;
% Evaluate the results.
fprintf('Rank of estimated matrix: %d \n', rank(Z));
fprintf('Relative Mean Square Error of estimated matrix with L1 regularization: %.2f%% \n',...
    (norm(Z(:) - actual_matrix(:))/norm(actual_matrix(:)))*100);

% Save the estimates and performance of algorithm.
imwrite(Z, strcat(filename, '/estimated_image_with_L1.png'));
formatSpec = 'Estimated image relative error with L1 regularization: %4.2f \r\n';
fprintf(fileID, formatSpec, (norm(Z(:) - actual_matrix(:))/norm(actual_matrix(:)))*100);
save(strcat(filename, '/all_variables/all_variables.mat'));

% Now the L1_norm is not included.
L1_norm = 0;
[X, U, Y, ~, ~] = admm(U_i, Y_i, observed_image, observed_image, Omega, L1_norm);

Z = X*U*Y;
% Evaluate the results.
fprintf('Rank of estimated matrix: %d \n', rank(Z));
fprintf('Relative Mean Square Error of estimated matrix without L1 regularization: %.2f%% \n',...
    (norm(Z(:) - actual_matrix(:))/norm(actual_matrix(:)))*100);

% Save the estimates and performance of algorithm.
imwrite(Z, strcat(filename, '/estimated_image_without_L1.png'));
formatSpec = 'Estimated image relative error without L1 regularization: %4.2f \r\n';
fprintf(fileID, formatSpec, (norm(Z(:) - actual_matrix(:))/norm(actual_matrix(:)))*100);
save(strcat(filename, '/all_variables/all_variables.mat'));
fclose(fileID);