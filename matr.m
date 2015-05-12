function Z = matr(v, n, r)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva
    
    % matr operator
    
    % v is a vector of size 1/2r(r+1)+(n-r) n
    Z = zeros(n, r);
    k = 1;
    for j = 1:r
        for i = j:n
            Z(i, j) = v(k);
            k = k + 1;
        end
    end
