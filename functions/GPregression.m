function[mean, covariance_matrix, K_X_X, inverse_time, mean_time, covariance_time,pre_allocation_time] = GPregression(X, y, Xstar, hyp)
    %defining sizes and constants given by input
    N_X = size(X, 1);
    N_Xstar = size(Xstar, 1);
    sigma_l = hyp(1); sigma_f = hyp(2);
    
    %First making a large kernel matrix of all the inputs (train+test)
    joint_kernel_matrix = kernel_matrix([X; Xstar], sigma_l, sigma_f);
    
    %sub-dividing the joint kernel matrix in its respective parts
    K_X_X = joint_kernel_matrix(1:N_X, 1:N_X);
    K_X_Xstar = joint_kernel_matrix(1:N_X,N_X+1:end);
    K_Xstar_X = K_X_Xstar.';
    K_Xstar_Xstar = joint_kernel_matrix(N_X+1:end, N_X+1:end);

    % Start specific timers with a handle to time calculation time of
    % various outputs.
    inv_timer = tic;
    K_X_X_inv = inv(K_X_X);
    inverse_time = toc(inv_timer);  % Stops the specific timer

    pre_timer = tic;
    pre = K_Xstar_X * K_X_X_inv;
    % clear("K_X_X_inv");
    % clear("K_Xstar_X");
    pre_allocation_time = toc(pre_timer);

    mean_timer = tic;
    mean = pre * y;
    % clear("y");
    mean_time = toc(mean_timer);

    cov_timer = tic;
    covariance_matrix = K_Xstar_Xstar - pre * K_X_Xstar;
    covariance_time = toc(cov_timer);
end
