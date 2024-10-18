function Kernel = kernel_matrix(X, sigma_l, sigma_f)
    % non-optimized version, Is correct (I think)
    N = size(X, 1);
    Kernel = zeros(N, N);
    for i = 1:N
        for j = i:N
            val = sigma_f * exp(-(norm(X(i,:) - X(j,:))^2)/(2*sigma_l));
            Kernel(i,j) = val;
            Kernel(j,i) = val;
        end
    end

    % % optimized version of the for loop
    % N = size(X, 1);
    % % permute X to [3,2,1] to make use of broadcasting when subtracting
    % X_permuted = permute(X, [3,2,1]);
    % % difference of each row with itself and every other row (in 3D)
    % Kernel = X_permuted - X;
    % % absolute value and summed over second dimension to obtain frobenius norm squared
    % Kernel = abs(Kernel);
    % Kernel = sum(Kernel, 2);
    % % rest of the squared exponential kernel function calculation
    % Kernel = squeeze(sigma_f .* exp(-0.5/sigma_l .* Kernel));
    
    %even new method
    %Kernel2 = sigma_f^2 .* exp(-0.5/sigma_l^2 .* pdist2(X, X).^2);

    %assert(sum(Kernel-Kernel2,'all') < 1e-9)

end

