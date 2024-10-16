function tensor = TT_reconstruct(tt_cores)
%TT_reconstruct(tt,n) 
%   Algorithm that reconstructs a N-th order tensor from the tensor train
%   decomposition
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
%   tensor (N-dimensional double):  N-th order tensor reconstructed from tt. 
%------ Implement your code below ------
N = length(tt_cores) - 1;  % Amount of cores in cell array. -1 because the last cell is the norm thiny.

% Initialization stuff
% if ndims(tt_cores{1}) == 2
%     Tr = mode_n_matricization(tt_cores{1}, 1);
% else
%     Tr = mode_n_matricization(tt_cores{1}, 2);
% end

Tr = mode_n_matricization(tt_cores{1},2);
sz_vec = zeros(1, N);
sz_vec(1) = size(Tr, 1);

for i = 2:N
    sz_vec(i) = size(tt_cores{i}, 2);
    sz(1) = size(tt_cores{i},1);
    sz(2) = size(tt_cores{i},2);
    sz(3) = size(tt_cores{i},3);
    Tr = Tr * mode_n_matricization(tt_cores{i}, 1);
    if i < N
        Tr = reshape(Tr, [], sz(3));  % Reshape only if not the last core
    end
end
% Reshape to get the multi-dimentional Tensor thingy
tensor = reshape(Tr, sz_vec);

end