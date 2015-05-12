function A = vec(Ainput)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

% vec operator
    [n, m] = size(Ainput);
    A = reshape(Ainput, n * m, 1);
end

