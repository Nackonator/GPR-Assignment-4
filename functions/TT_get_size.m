function sz = TT_get_size(tt)
%TT_get_size (tt,n) 
%   Algorithm that returns the size of a tt
%      
%INPUT:
%   tt (cell array with N+1 cells):  tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0.
%OUTPUT:
%   sz (N x 3 double):              size of TT decomposition, see example
%                                   figure 
%------ Implement your code below ------
N = numel(tt)-1;
sz = zeros(N,3);

%I noticed that some shapes were 3D and other were 2D, so I just forced
%each dimention to be assigned a size (even if it's 1).
for i = 1:N
    sz(i,1) = size(tt{i},1);
    sz(i,2) = size(tt{i},2);
    sz(i,3) = size(tt{i},3);
end
end