function [X,y,Xstar] = getdata(Nd,hyp)
    
    coordX1    = linspace(-1,1,Nd)';
    coordX2    = linspace(-1,1,Nd)';
    coordX3    = linspace(-1,1,Nd)';
    [X1,X2,X3] = meshgrid(coordX1,coordX2,coordX3);
    X          = [X1(:),X2(:),X3(:)];

    K = zeros(size(X,1),size(X,1));
    for i = 1:size(X,1)
        for j = 1:size(X,1)
            exparg = norm(X(i,:)-X(j,:));
            K(i,j) = hyp(2) * exp(-(exparg^2/(2*hyp(1))));
        end
    end
    
    f          = chol(K+1e-5*eye(size(K)))*randn(Nd*Nd*Nd,1);
    y          = f;
    [X1,X2,X3] = meshgrid(linspace(-1.1,1.1,4*Nd),linspace(0,0,1),linspace(-1.1,1.1,4*Nd));
    Xstar      = [X1(:),X2(:),X3(:)];
    
end

