function out = kron_mvprod(K,b)
%K should be a cell array with all the K1, K2 etc
%b should be the vector you want to multiply the K(X,X)^-1 with.

D = length(K);
N = length(b);

x = b;

for d = 1:D
    Gd = length(K{d});
    X = reshape(x,Gd,N/Gd);
    Z = K{d}\X;
    Z = Z';
    x = vec(Z);
end
out = x;
end