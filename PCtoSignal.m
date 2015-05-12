function [outputSignal, outputCovariances] = PCtoSignal(inputSignalPC, inputCovPC, V)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % reproduces the outputSignal out of the principal components provided
    % in inputSignalPC
    
    % V is the matrix containing the eigenvalues of the cov matrix in the columns
    if ~isempty(inputSignalPC)
        outputSignal = inv(V) * inputSignalPC;
    else
        outputSignal = [];
    end
    
    if ~isempty(inputCovPC)
        [T, N, ~] = size(inputCovPC);
        outputCovariances = zeros(T, N, 1);
        invV = inv(V);
        invVT = invV';
        for i = 1:T
            htemp = math(reshape(inputCovPC(i, :, 1), N, 1));
            outputCovariances(i, :, 1) = vech(invV * htemp * invVT)';
        end 
    else
        outputCovariances = [];
    end
    