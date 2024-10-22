clc
TT_Kkron{1} = inv(K3);
TT_Kkron{2} = inv(K2);
TT_Kkron{3} = inv(K1);
s = randn(size(y1));


Kkron_inv = kron(inv(K3),kron(inv(K2),inv(K1)));

test1 = Kkron_inv*y1;
test2 = kron_mvprod(TT_Kkron,y1);
norm(test1-test2,'fro')


test3 = Kkron_inv * [y1 s];
test4 = kron_mmprod(TT_Kkron,[y1 s]);
norm(test3-test4,'fro')