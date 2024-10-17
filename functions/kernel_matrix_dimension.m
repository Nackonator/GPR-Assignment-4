%coding equation 5.4
function K = kernel_matrix_dimension(x, sigma_l, sigma_f)
K = zeros(length(x));

for i = 1:length(x)
    for j = 1:length(x)
        K(i,j) = sigma_f^2*exp(-0.5/sigma_l^2 * (x(i)-x(j))^2);
    end
end