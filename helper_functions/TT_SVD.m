function [tt_cores,rel_error] = TT_SVD(tensor, epsilon)
%TT_SVD(tensor,epsilon) 
%   Algorithm that decomposes a given N-th order tensor into tensor train format
%INPUT:
%   tensor (N-dimensional double):  N-th order tensor
%   epsilon (double):               error-bound for approximation error of TT decomposition
%OUTPUT:
%   tt_cores (cell array with N+1 cells): tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0.
%   error (double):                 actual approximation error of TT
%                                   decomposition
n = ndims(tensor);
% every cell(1,n) contains one tt-core, which is a r_i x n_i x r_i+1 double
% do not change the datastructure! Otherwise, the testfunction is not
% guaranteed to work.
% The cell(1,n+1) contains the location of the core at the moment. 
% If the norm is stored in the 10th core, the value in the cell(1,d+1) should be 10. 
% If unsure, have a look at the variable tt_dog
tt_cores = cell(1,n+1);

%------ Implement your code below ------
sz = size(tensor);
Tr = mode_n_matricization(tensor,1);
R_list = zeros(1,n);
rel_error = 0;
frb = frob_norm(tensor);

for i = 1:n-1
    R = 0;          %Reset rank cut-off for new SVD decomp.
    [U,S,V] = svd(Tr,"econ");
    err = frob_norm(S);

    % Finding truncation/R_i variable which allows the error bound thingy
    while err^2 > (epsilon^2*frb^2)/(n-1)
        R = min( [R+1 , prod(sz(1:i)) , prod(sz(i+1:end))] ); %start with R = 0+1 = 1 seeing as you need atleast R=1 for the SVD to even work.
        S2 = S( R+1:end , R+1:end );
        err = frob_norm(S2);
    end

    % Final check to see if R isn't too large
    if R > 1
        S2_prev = S(R:end, R:end);
        err_prev = frob_norm(S2_prev);
        if err_prev^2 <= (epsilon^2 * frb^2) / (n-1)
            R = R - 1;
            err = err_prev; % Update to the correct error
        end
    end

    R_list(i) = R;
    rel_error = rel_error + err^2;

    %Defining truncated SVD matrices
    U1 = U(:,1:R);
    S1 = S(1:R,1:R);
    V1 = V(:,1:R);

    %reshaping to get tt_cores and next Tr matrix
    Tr = S1*V1';
    Tr = reshape(Tr,R_list(i)*sz(i+1),[]);
    
    %calculating the tt_cores.
    if i == 1
        tt_cores{i} = reshape(U1,1,sz(i),R_list(i));
    else
        tt_cores{i} = reshape(U1,R_list(i-1),sz(i),R_list(i));
    end

end
tt_cores{n} = S1*V1';
tt_cores{n+1} = n;%norm(mode_n_matricization(tt_cores{n},2));

rel_error = sqrt(rel_error)/frb;
end