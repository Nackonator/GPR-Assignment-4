function tt = site_n(tt,n)
%site_n(tt,n) 
%   Algorithm that brings tensor train into site_n mixed canoncical form. 
%   decomposition 
%      
%INPUT:
%   tt (cell array with N+1 cells): tensor train (tt) decomposition of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell. If the tt is not in site_n
%                                   mixed canonical form, the N+1 cell
%                                   contains a 0. 
%   n (int):                        n with 1<=n<=N, location where norm is 
%                                   supposed to be moved to.
%OUTPUT:
%   tt (cell array with N+1 cells): tensor train (tt) decomposition with 
%                                   norm in positin n. It is a tt of a 
%                                   N-th order tensor, where the tt-cores 
%                                   are stored in the first N cells and the 
%                                   location of the norm-core is stored in the 
%                                   N+1 cell.
%------ Implement your code below ------

N = length(tt)-1;

%Left to Right part
for i = 1:n-1
    %Saving size for reshaping
    sz1 = size(tt{i});

    %Reshaping Tensor in order to apply QR decomposition
    if ismatrix(tt{i})
        Tr = tt{i};
    % elseif i==1
    %     Tr = mode_n_matricization(tt{i},2);
    else
        Tr = mode_n_matricization(tt{i},3);
    end
    Tr = Tr';
    [Q,R] = qr(Tr,'econ');
    
    % S_Tr = size(Tr)
    % Sr = size(R)
    % Sq = size(Q)
    % Stt2 = size(tt{i+1})

    tt{i}   = reshape(Q,sz1(1),sz1(2),[]);
    tt{i+1} = mode_n_product(tt{i+1},R,1);

end

%Right to left part
for j = 0:N-n-1
    %j
    %Saving size for reshaping
    sz2 = size(tt{N-j});
        
    %Doing LU decomp but the with the qr() command
    if ismatrix(tt{N-j})
        Tr = tt{N-j}';
    else
        Tr = mode_n_matricization(tt{N-j},1);
        Tr = Tr';
        %Tr = reshape(tt{N-j});
    end

    [Q,R] = qr(Tr,'econ');
    L = R';
    Q = Q';

    % Size_Tr = size(Tr)
    % Size_L = size(L)
    % Size_Q = size(Q)
    % Size_tt2 = size(tt{N-j-1})
    % N-j-1

    %re-defining tensors and stuff
    if j == 0
        tt{N-j} = reshape(Q,[],sz2(2));
    else
        tt{N-j} = reshape(Q,[],sz2(2),sz2(3));
    end
    tt{N-j-1} = mode_n_product(tt{N-j-1},L',3); %from left so now mode-3-product.
end
tt_n = n;%tt{n};
tt{N+1} = frob_norm(tt_n);

end