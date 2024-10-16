N_X = size(X, 1);
N_Xstar = size(Xstar, 1);

% fit GP model
%GPmodel = fitrgp(X, y, 'KernelFunction', 'squaredexponential', 'KernelParameters', sqrt(hyp));
% Get fitted Kernel param values
%hyp2 = GPmodel.KernelInformation.KernelParameters;

sigma_l = hyp(1); sigma_f = hyp(2);
disp('hyper parameters defined')
%predictions = predict(GPmodel, Xstar);

joint_kernel_matrix = kernel_matrix([X; Xstar], sigma_l, sigma_f);
disp('joint kernel matrix calculated')

% getting different parts of the joint covariance matrix
K_X_X = joint_kernel_matrix(1:N_X, 1:N_X);
K_X_Xstar = joint_kernel_matrix(1:N_X,N_X+1:end);
K_Xstar_X = K_X_Xstar.';
K_Xstar_Xstar = joint_kernel_matrix(N_X+1:end, N_X+1:end);
disp('k_x matrices made')

mean = zeros(size(K_Xstar_X, 1), size(y, 2));
covariance_matrix = zeros(size(K_Xstar_Xstar));
disp('data pre-allocated')

%DEBUGGING STUFF, IGNORE OR REMOVE LATER
% size_K_X_X          = size(K_X_X)
% size_K_X_Xstar      = size(K_X_Xstar)
% size_K_Xstar_X      = size(K_Xstar_X)
% size_K_Xstar_Xstar  = size(K_Xstar_Xstar)

% calculate mean vector and covariance matrix based on eq 2.19 from the book
%--+++====BACKSLASH OPERATOR ISNIET MOGELIJK OMDAT DE NUMBER OF ROWS NIET HETZELFDE ZIJN.====+++--
tic
K_X_X_inv = inv(K_X_X); %pre-calculate the inverse seeing as it is needed in both the mean and covariance. This is also a very intensive inverse to calculate (big matrix).
inverse_time = toc

disp('inverse done')

tic
pre = K_Xstar_X * K_X_X_inv;
pre_allocation_time = toc

disp('pre-allocation of re-used matrix in mean and covariance calculation done')

tic
mean = pre * y;
mean_time = toc
disp('mean done')

tic
covariance_matrix = K_Xstar_Xstar - pre*K_X_Xstar;
covariance_time = toc

disp('mean done')