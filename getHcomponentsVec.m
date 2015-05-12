function [ HA, HB] = getHcomponentsVec(h, indArray)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % get the blocks of the Hessian for A and B matrices
    % in vec form
    HA = (h(1:indArray(1)));

    HB = (h(indArray(1) + 1:indArray(2))); 
        
end

