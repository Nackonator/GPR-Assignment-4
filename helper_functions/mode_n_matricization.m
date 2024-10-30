function X = mode_n_matricization(X,n)
% MODE_K_MATRICIZATION takes a tensor X as input and a mode n such
% that n<ndims(X) and returns the mode-n matricization of X.
% INPUTS tensor X, mode n.
% OUPUT mode-n matricization of X.

% First "measuring/determining" size and dimension of Tensor
X_size = size(X);           %Get the I values (e.g. how big each dimension is)
Xn_dim = X_size(n);         %Check the size\length of the n-th dimension
X_dim = length(X_size);     %Determine the Dimension of the tensor (e.g. is it 3D, 4D etc.)
X_ind = 1:X_dim;

% Changing "order" of tensor
X_newind = [X_ind(n),X_ind(1:n-1),X_ind(n+1:end)]; %reshuffle indices

% Reshaping the "shuffled" tensor
X_perm = permute(X,X_newind);    %Permutate them bitches (respectfully)
X = reshape(X_perm,Xn_dim,[]);   %Reshape it into a Matrix
end 