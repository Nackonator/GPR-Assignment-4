function out = kron_mmprod(K,B)
% K should be a cell array with all the K1, K2, etc.
% B is the matrix you want to multiply with the Kronecker product.
% out is the result of kron(K1, ..., Kd) * B

D = length(K);      % Number of matrices
[N,M] = size(B);   % Dimensions of the input matrix B
x = B;              % Initialize x with input matrix B

for d = D:-1:1      % Loop backwards from D to 1
    Gd = length(K{d});
    for m = 1:M
        X = reshape(x(:,m),Gd,N/Gd);
        Z = K{d}\X;
        Z = Z';
        x(:,m) = vec(Z);
    end
end
out = x;
end
