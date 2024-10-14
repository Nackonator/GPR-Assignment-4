function[GPmodel, predictions, mean, covariance_matrix, Kernel_train_data] = GPregression(X, y, Xstar, hyp)
    N_X = size(X, 1);
    N_Xstar = size(X, 1);

    % fit GP model
    GPmodel = fitrgp(X, y, 'KernelFunction', 'squaredexponential', 'KernelParameters', sqrt(hyp));
    % Get fitted Kernel param values
    [sigma_l, sigma_f] = GPmodel.KernelInformation.KernelParameters;

    % get return values
    Kernel_train_data = kernel_matrix(X, sigma_l, sigma_f);
    predictions = predict(GPmodel, Xstar);

    joint_kernel_matrix = kernel_matrix([X; X_star], sigma_l, sigma_f);

    % getting different parts of the joint covariance matrix
    K_X_X = joint_kernel_matrix(1:N_X, 1:N_X);
    K_X_Xstar = joint_kernel_matrix(N_X + 1:N_X + N_Xstar);
    K_Xstar_X = K_X_Xstar.';
    K_Xstar_Xstar = joint_kernel_matrix(N_X + 1:N_X+N_Xstar, N_X + 1:N_X + N_Xstar);

    % calculate mean vector and covariance matrix based on eq 2.19 from the book
    mean = K_X_X \ K_Xstar_X .* predictions;
    covariance_matrix = K_Xstar_Xstar - K_X_X \ K_Xstar_X * K_X_Xstar;
end