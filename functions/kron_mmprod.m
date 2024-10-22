function out = kron_mmprod(K,B)
% K should be a cell array with all the K1, K2, etc.
% B is the matrix you want to multiply with the Kronecker product.
% out is the result of kron(K1, ..., Kd) * B

D = length(K);      % Number of matrices
[M, N] = size(B);   % Dimensions of the input matrix B
x = B;              % Initialize x with input matrix B

for d = D:-1:1      % Loop backwards from D to 1
    Gd = size(K{d}, 1);         % Get the size of the current matrix
    numCols = N / Gd;           % Number of columns after reshaping
    
    % Initialize Z to store the result for each column
    Z = zeros(Gd, numCols * M / Gd); 
    
    for i = 1:M  % Loop over each column of B
        X = reshape(x(:, i), Gd, numCols);  % Reshape the column of x
        Z_col = K{d} * X;                      % Matrix multiplication
        Z(:, i) = Z_col(:);                 % Store the result as a vector
    end
    
    x = Z';  % Transpose and assign the result to x
end

out = x;   % The final result is stored in x
end
