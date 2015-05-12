function A = vechstar(v)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % vech* operator
    
    A = 0.5 * (math(v) + diag(diag(math(v))));