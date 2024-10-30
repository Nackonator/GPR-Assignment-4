function Sum = TT_add(A, B)
%TT_add(A, B) 
%   Adds two tensors in tensor train (TT) format.
%
%INPUT:
%   A, B: Cell arrays representing tensor train decompositions.
%OUTPUT:
%   Sum: A cell array representing the sum of the tensors in TT format.

N = length(A) - 1;  % Number of cores (ignoring the norm-core)
Sum = cell(1, N + 1);  % Initialize the output tensor train

for i = 1:N
    sz_A = size(A{i});
    sz_B = size(B{i});
    
    %test to see if 2nd dimention (I_n dimentions) is correct/equal
    assert(sz_A(2) == sz_B(2),'DIM MISMATCH')
    
    if i == 1
        % First core: horizontal concatenation
        Sum{i} = [reshape(A{i},sz_A(2),sz_A(3)) reshape(B{i},sz_B(2),sz_B(3))];
        sz_Sum = size(Sum{i});
        Sum{i} = reshape(Sum{i},[1,sz_Sum]);
        
    elseif i == N
        % Last core: vertical concatenation
        Sum{i} = [A{i}; B{i}];
        
    else

        % Method 1
        % % Allocate with zeros
        % Sum{i} = zeros(sz_A(1) + sz_B(1), sz_A(2), sz_A(3) + sz_B(3));
        % 
        % % Fill the block diagonal
        % Sum{i}(1:sz_A(1), :, 1:sz_A(3)) = A{i};
        % Sum{i}(sz_A(1) + 1:end, :, sz_A(3) + 1:end) = B{i};

        %Method 2, as in lecture
        Sum{i} = zeros(sz_A(1) + sz_B(1), sz_A(3)+sz_B(3), sz_A(2)); %allocating a zero tensor
        T_A = permute(A{i},[1,3,2]);
        T_B = permute(B{i},[1,3,2]);

        Sum{i}(1:sz_A(1),1:sz_A(3),:) = T_A;
        Sum{i}(sz_A(1)+1:end,sz_A(3)+1:end,:) = T_B;

        Sum{i} = permute(Sum{i},[1,3,2]);
    end
end

% Handle the norm-core (last element in the cell array)
Sum{end} = 10;  % Assuming that neither tensor is in the mixed canonical form

end
