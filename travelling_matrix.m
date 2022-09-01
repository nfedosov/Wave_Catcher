function [M] = travelling_matrix(X,dX)


    koef = norm(X,'fro');
    X = X./koef;
    dX = dX./koef;
    Nch = size(X, 2);
    M = ones(Nch, Nch)/Nch;

    Niter = 100000;
    rate = 0.001;
    lam = 1;
    d2 = 0;
    for i = 1:Niter
        d1 = real(X'*(X*M-dX));
        %d2 = real(lam*(sign(M.*M').*M'));
        %for j = 1:Nch
        %    d2(j,j) = 0;
        %end
        
        M = M - rate*(d1+d2);
        %M = M.*eye(Nch);
        M = (M+M')/2;
        %M = (M-M')/2;
    end
a = 1;
end