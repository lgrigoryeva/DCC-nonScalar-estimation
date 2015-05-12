function A = Projection_nS(B)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % Computes the n-symmetric projection of a n^2xn^2 matrix
    n = round(sqrt(size(B, 1)));
    for k = 1:n
       for l = 1:n
            temp = B(n * (k - 1) + 1:k * n, n * (l - 1) + 1:l * n);
            temp = .5 * (temp + temp');
            B(n * (k - 1) + 1:k * n, n * (l - 1) + 1:l * n) = temp;
       end
    end
    A = .5 * (B + B');