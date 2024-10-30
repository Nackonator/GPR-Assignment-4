function Kernel = GPR_GRID_kernel2(x, sigma_l, sigma_f)
   Kernel =  x - x';
   Kernel = sigma_f .* exp(-0.5/sigma_l .* Kernel.^2);
end

