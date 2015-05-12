function res = vech(A)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % vech function that stacks the lower triangular portion of a n*n matrix as
    % a n*(n+1)/2 vector. Its inverse is the math function
    
    [~, n] = size(A);
    res = zeros(n * (n + 1) / 2, 1);

    tempbound = 0;
    for i = 1:n
       res(tempbound + 1:tempbound + n - i + 1, 1) = A(i:n, i); 
       tempbound = tempbound + (n - i) + 1;
    end