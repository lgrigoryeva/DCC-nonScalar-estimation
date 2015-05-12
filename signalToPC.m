function [outputData, newV] = signalToPC(inputData, numPrincipalComp, V)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva
    
    % provides the PCA for the inputData
    centeredInputData = inputData;
    if ~isempty(V)
        outputData = V * centeredInputData;
        newV = V;
    else
        covMatrix = cov(centeredInputData');

        [V, D] = eig(covMatrix);
        eigenvals = diag(D);
        [newD, indexes] = sort(eigenvals, 'descend');
        newV = V(:, indexes(1:numPrincipalComp));
        outputData = diag(sqrt(1./(newD))) * newV' * centeredInputData;
        newV = diag(sqrt(1./(newD))) * newV';
    end
    
