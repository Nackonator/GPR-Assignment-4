function Kernel = axis_aligned_kernel(x, sigma_l, sigma_f)
   Kernel =  x - x';
   Kernel = sigma_f^2 .* exp(-0.5/sigma_l^2 .* Kernel);
end

