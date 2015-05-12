function dimParam = getModelNumParameters(modelTypeName, n, modelFamily)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

% computes the number of the intrinsic parameters depending on the model
% specification

    if strcmp(modelTypeName,'hadamard')
        dimParam = n * (n + 1);
    elseif strcmp(modelTypeName,'scalar')
        dimParam = 2;
    elseif strcmp(modelTypeName,'rankdef1')
        r = 1;
        dimParam = 2 * (n * r - 0.5 * r * (r - 1));
    elseif strcmp(modelTypeName,'rankdef2')
        r = 2;
        dimParam = 2 * (n * r - 0.5 * r * (r - 1));
    elseif strcmp(modelTypeName, 'almon')
        dimParam = 6;
    end
    % by default modelFamily = 'dcc'
    if nargin < 3
        dimParam = dimParam + 3 * n;
    end
end

