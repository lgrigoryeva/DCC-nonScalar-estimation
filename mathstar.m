function v = mathstar(A)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % Adjoint of the math operator
    v = 2 * vech(A - .5 * diag(diag(A)));