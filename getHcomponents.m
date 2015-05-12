function [HA, HB] = getHcomponents(h, indArray)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % get the blocks of the Hessian for A and B matrices
    HA = (h(1:indArray(1)));

    HB = mat(h(indArray(1) + 1:indArray(2))); 
        
end

