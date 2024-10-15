function Kernel = kernel_matrix_unopt(X, sigma_l, sigma_f)
    % non-optimized version
    N = size(X, 1);
    Kernel = zeros(N, N);
    size(Kernel)
    for i = 1:N
        for j = i:N
            val = sigma_f^2 * exp(-(norm(X(i,:) - X(j,:))^2)/(2*sigma_l^2));
            Kernel(i,j) = val;
            Kernel(j,i) = val;
        end
    end
    size(Kernel)
end
