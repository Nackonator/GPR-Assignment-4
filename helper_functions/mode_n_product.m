function Z = mode_n_product(X,Y,n)
    % MODE_N_PRODUCT takes tensor X and compatible matrix Y and performs mode-n product between X and Y.
    % INPUT tensor X, matrix Y.
    % OUTPUT tensor Z.

    % Getting data from input, e.g. sizes, dimensionality etc. and also apllying the mode-n matricization thingy
    X_modn = mode_n_matricization(X,n);
    size_X = size(X);
    size_Y = size(Y);
    X_dim = length(size_X);
    X_ind = 1:X_dim;
    
    % Repeating some stuff I did in the mode_n_matricization
    X_newind = [X_ind(n),X_ind(1:n-1),X_ind(n+1:end)]; %reshuffle indices
    X_perm = permute(X,X_newind); %Permute the X tensor in order to get the correct dimension sizes
    size_X_perm = size(X_perm);   %Getting the correct dimension sizes

    % Calculating the matrix product of the n-mode matricization thingy
    Z_2d = Y*X_modn;

    % Folding it back into the original "dimensionality" of the tensor (e.g. if the tensor was a 3D tensor, the Z_2d will be "folded" back into a 3D shape (with changed sizes)
    size_X_perm(1) = size_Y(1);      %Need to "change" shape of tensor back to its original dimension. To do that the size of the dimensions needs to be determined.
    Zt = reshape(Z_2d,size_X_perm);  %Changing it back to its dimension
    
    % Permuting Z back e.g. 'undoing the permutation nececarry for the
    % n-mode matricization
    X_ind2 = 1:X_dim;
    X_ind2(1) = [];
    X_ind2 = [X_ind(2:n),1,X_ind(n+1:end)];
    Z = permute(Zt,X_ind2);
end