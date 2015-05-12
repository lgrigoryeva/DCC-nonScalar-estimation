function v = vecr(A)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva
    
    % vecr operator
    [n, r] = size(A);
    v = zeros(n * r - 0.5 * n * (n - 1), 1);
    
    tempbound = 1;
    
    for i = 1:r
       v(tempbound:tempbound + n - i, 1) = A(i:n, i); 
       tempbound = tempbound + (n - i + 1);
    end    

end

