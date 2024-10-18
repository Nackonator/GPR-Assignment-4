%coding equation 5.4
function K = kernel_matrix_dimension(x, sigma_l, sigma_f)
K = zeros(length(x));

for i = 1:length(x)
    for j = i:length(x)
        K(i,j) = sigma_f*exp(-(x(i)-x(j))^2/(2*sigma_l));
        K(j,i) = sigma_f*exp(-(x(i)-x(j))^2/(2*sigma_l));
    end
end