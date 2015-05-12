function Z = mat(v)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % Inverse of the vec operator
    % Given a (vertical) n^2 dimensional vector v, z is the corresponding
    % nxn  matrix obtained by filling z by columns

    n = round(sqrt(length(v)));
    Z = reshape(v, n, n);

