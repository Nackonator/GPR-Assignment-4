function tt_cores = TT_round(tt,epsilon)
%TT_round (tt,n) 
%   Algorithm that rounds a TT algorithm up to a defined approximation
%   error
%      
%INPUT:
%   tt (cell array with N+1 cells): tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0.
%   epsilon (double):               error-bound for approximation error of TT decomposition
%OUTPUT:
%   tt_cores (cell array with N+1 cells): tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0.
%------ Implement your code below ------
N = length(tt) - 1;
tt = hidden_site_n(tt,N); %I'm making sure my tensor train really has its norm in the very last tensor core. I realised halfway that it's maybe more efficient to just set this to site-1 and then do the svd from left to right but this code works so I'm not chaning it.
tt_cores = cell(1,N+1);

n = tt{end};

Tr = tt{N};
S_Tr = size(Tr);
frb = frobnorm_tt(tt);
R_list = zeros(1,N);
for i = 0:N-2
    sz = size(tt{N-i});
    sz2 = size(tt{N-i-1});
    R = 0;
    [U,S,V] = svd(Tr,'econ');
    err = frob_norm(S);

    while err^2 > (epsilon^2*frb^2)/(N-1)
        R = R+1; %Not adding size checks cuz idk how to do it when it's in Tensor train format
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
    R_list(i+1) = R;
    U1 = U(:,1:R);
    S1 = S(1:R,1:R);
    V1 = V(:,1:R);
    
    if ismatrix(tt{N-i})
        tt_cores{N-i} = V1';
    else
        %size(V1')
        tt_cores{N-i} = reshape(V1',R_list(i+1),sz(2),R_list(i));
    end
    tt_cores{N-i-1} = mode_n_product(tt{N-i-1},(U1*S1)',3);
    Tr = reshape(tt_cores{N-i-1},sz2(1),[]);
end

tt_cores{1} = mode_n_product(tt{1},(U1*S1)',3);
tt_cores = hidden_site_n(tt_cores,10);


tt_cores{end} = 10;
end