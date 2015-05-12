function indexes = reorderToPC(inputData, numPrincipalComp)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % finds the order in inputData w.r.t. the descending magnitude of the
    % projections on the first principal component
    covMatrix = cov(inputData');

    [V, D] = eig(covMatrix);
    eigenvals = diag(D);
    [~, indexes] = sort(eigenvals, 'descend');
    newV = V(:, indexes(1:numPrincipalComp));
    outputData = newV' * inputData;

    importArray = abs(inputData * outputData(1, :)');
    [~, indexes] = sort(importArray, 'descend');

    
