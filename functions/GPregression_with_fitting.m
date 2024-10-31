function[GPmodel, predictions, mean, covariance_matrix, Kernel_train_data] = GPregression_with_fitting(X, y, Xstar, hyp)
    N_X = size(X, 1);
    N_Xstar = size(Xstar, 1);

    % fit GP model
    GPmodel = fitrgp(X, y, 'KernelFunction', 'squaredexponential', 'KernelParameters', sqrt(hyp));

    % Get fitted Kernel param values
    hyp2 = GPmodel.KernelInformation.KernelParameters;
    sigma_l = hyp2(1); sigma_f = hyp2(2);

    % get return values
    Kernel_train_data = kernel_matrix(X, sigma_l, sigma_f);
    predictions = predict(GPmodel, Xstar);

    joint_kernel_matrix = kernel_matrix([X; Xstar], sigma_l, sigma_f);

    % getting different parts of the joint covariance matrix
    K_X_X = joint_kernel_matrix(1:N_X, 1:N_X);
    K_X_Xstar = joint_kernel_matrix(1:N_X,N_X+1:end);
    K_Xstar_X = K_X_Xstar.';
    K_Xstar_Xstar = joint_kernel_matrix(N_X+1:end, N_X+1:end);

    %DEBUGGING STUFF, IGNORE OR REMOVE LATER
    % size_K_X_X          = size(K_X_X)
    % size_K_X_Xstar      = size(K_X_Xstar)
    % size_K_Xstar_X      = size(K_Xstar_X)
    % size_K_Xstar_Xstar  = size(K_Xstar_Xstar)

    % calculate mean vector and covariance matrix based on eq 2.19 from the book
    K_X_X_inv = inv(K_X_X); %pre-calculate the inverse seeing as it is needed in both the mean and covariance. Done seperately in order to time it.
    mean = K_Xstar_X * K_X_X_inv * y;
    covariance_matrix = K_Xstar_Xstar - K_Xstar_X*K_X_X_inv*K_X_Xstar;
end