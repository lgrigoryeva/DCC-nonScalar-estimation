function res = math(v)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % Inverse of the vech operator
    % Given a (vertical) 1/2(n+1)n dimensional vector v, z is the corresponding
    % nxn symmetric matrix obtained by filling z by rows
    
    u = v';
    n = round(-.5 + .5 * sqrt(1 + 8 * length(v)));
    
    res = zeros(n);
    
    j = 1;
    for i = 1:n    
       res(i, i:n) = u(j:j + n - i); 
       j = j + n - i + 1;
    end
    
    res = res + res' - diag(diag(res));