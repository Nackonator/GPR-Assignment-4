function[pred_mean, pred_cov] = GPregression(X,y,Xstar,hyp)

    N = dim(X,1);

    GPmodel = fitrgp(X,y,'KernelFunction','squaredexponential','KernelParameters',hyp); %Calculate the GP model
    [l,sign] = GPmodel.KernelInformation.KernelParameters; %!!! Check if square is nececary !!!

    K_matrix = zeros(N,N); %allocate for speed

    for i = 1:N
        for j =i:N
            val = sign^2*exp(-(norm(X(i,:)-X(j,:),'fro')^2)/(2*l^2));
            K_matrix(i,j) = val;
            K_matrix(j,i) = val;
        end
    end



end