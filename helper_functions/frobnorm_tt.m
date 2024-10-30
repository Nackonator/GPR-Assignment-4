function tt_norm = frobnorm_tt(tt_site_n)
%frobnorm_tt(tt,n) 
%   Algorithm that computes the frobenius norm of a TT, which is in site-n mixed canonical form 
%      
%INPUT:
%   tt (cell array with N+1 cells): tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. The TT is in site-n mixed
%                                   canonical form
%OUTPUT:
%   tt_norm (double):               frobenius norm of corresponsing tensor
%------ Implement your code below ------
N = length(tt_site_n)-1;

if ismatrix(tt_site_n{1})
    T = tt_site_n{1};
else
    T = mode_n_matricization(tt_site_n{1},2);
end

Mat = T'*T;

for i = 2:N-1
    T_up  = mode_n_product(tt_site_n{i},Mat',1);
    T_low = tt_site_n{i};

    Mat = mode_n_matricization(T_up,3) * mode_n_matricization(T_low,3)';
end
Mat = Mat*tt_site_n{N};

Mat = reshape(Mat,1,[]);
T2 = reshape(tt_site_n{N},[],1);

tt_norm = sqrt(Mat*T2);