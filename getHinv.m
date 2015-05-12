function Hinv = getHinv(H)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva
    
    % computes the inverse with SVD
    [U, S, V] = svd(H);
    Hinv = U * diag(1./(abs(diag(S)))) * V';
end

