function [Coefficients, LLF, Innovations, Sigmas ,Epsilons, LogL] = garch11(z, Coefficients)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

% Fits a GARCH(1,1) model with Gaussian innovations to the time series contained in the rows of z. 
% Preliminary estimation is carried out by fitting an ARMA(1,1) model to the squared time series.
% Coefficients is the array of structures with the fields that correspond
% to the ones we want to fit the sample to
[n, T] = size(z);

LLF = zeros(n, 1);
Innovations = zeros(T, n);
Sigmas = zeros(T, n);
Epsilons = zeros(T, n);

z = z';

if nargin < 2
    Coefficients = struct('coef', []);
    opts1 = struct(optimoptions('fmincon','Display','off', 'TolFun',1e-10));   
    garchModel = garch('GARCHLags',2,'ARCHLags',1, 'Distribution', 'Gaussian');
    %garchModel = garch('GARCHLags',1,'ARCHLags',1, 'Distribution', 't');
    for i = 1:n
        Coefficients(i) = struct('coef', []);
        [Coefficients(i).coef, ~, LogL, info] = estimate(garchModel, z(:, i), 'options', opts1);
    end
end

for i = 1:n
    [Sigmas(:,i), LLF(i)] = infer(Coefficients(i).coef, z(:, i));
    Sigmas(:,i) = sqrt(Sigmas(:,i));
    Innovations(:,i) = z(:,i) - Coefficients(i).coef.Offset;
    Epsilons(:,i) = Innovations(:,i)./Sigmas(:,i);
end

Epsilons = Epsilons';
Innovations = Innovations';
Sigmas = Sigmas';
    