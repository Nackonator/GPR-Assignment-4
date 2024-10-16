function z = frob_inner(X,Y)
    % INNER takes as input tensor X and tensor Y of equal sizes and of any order and returns 
    % their Frobenius inner product.
    % INPUT tensor X, Y.
    % OUTPUT scalar z.
    
    X_vec = reshape(X,1,[]);
    Y_vec = reshape(Y,[],1);

    z = conj(X_vec)*Y_vec;
end