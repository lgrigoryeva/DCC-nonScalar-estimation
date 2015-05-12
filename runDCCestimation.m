% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva
clear all
% turn the warnings off
warning('off','all');
% load the data
load('exampleData.mat');
% the file contains the series Returns_all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimension of the asset log-returns series
n = 5;

epsProjection_S = 1e-4;
% set the ascending sequence of the assets' numbers
combinCurrent = 1:n;

% choose the sample length reserved for the estimation T_est
[~, T_est] = size(Returns_all);
% initialize the sample structure
Sample = struct('series', Returns_all(combinCurrent, 1:T_est), ...
    'H0', cov(Returns_all(combinCurrent, 1:T_est)'), ...
    'z0', Returns_all(combinCurrent, 1), 'h', []);

ModelOutGarchFit = struct('GARCH_Coeffs', []);

[ModelOutGarchFit.GARCH_Coeffs, ~, ~, Sample.sigmas, Sample.Epsilons] = garch11(Sample.series);
S0 = cov(Sample.Epsilons');

% set up the flag whether to switch in the plotting
ifShowPlots = 1;

% set the estimation structure:
ModelOutLHD = struct('outStruct_scalar', [], 'outStruct_hadamard', [], ...
    'outStruct_rankdef1', [], 'outStruct_rankdef2', [], ...
    'outStruct_almon', [], 'effortStat_scalar', [], ...
    'effortStat_hadamard', [], 'effortStat_rankdef1', [], ...
    'effortStat_rankdef2', [], 'effortStat_almon', [], 'combin', combinCurrent);

% dimension of the parameter subspace for the Hadamard case
N = .5 * n * (n + 1);

% setting up the upper bound for the entries
K = 50;

%%%%%%%%%%%%%%%%%%%
% DCC scalar %
%%%%%%%%%%%%%%%%%%%

% estimation of the scalar model
% the scalar model can be estimated without using the Bregman divergencies
% used only as a benchmark

% setting the initial point
a_initialValue = 0.2;
b_initialValue = 0.7;

ModelDCC_scalar = struct('a', a_initialValue, 'b', b_initialValue, 'S', S0, ...
    'Q', S0, 'modelType', 'scalar', 'n', n, 'N', N);

% check whether the parameter constraints are satisfied for the chosen
% initial values
constraintsStruct_scalar = getConstraintVars(ModelDCC_scalar);

arrayOfConstraintsVals_scalar = getConstraintsCheckVals(constraintsStruct_scalar, ModelDCC_scalar);

% replace the initial values by the random ones if the constraints are not
% satisfied
% the initial random values are generated until the ones that satisfy the
% constraints are found

while (min(arrayOfConstraintsVals_scalar) <= epsProjection_S)
    
    ModelDCC_scalar.a = randn(1, 1);
    ModelDCC_scalar.b = randn(1, 1);
    
    constraintsStruct_scalar = getConstraintVars(ModelDCC_scalar);
    arrayOfConstraintsVals_scalar = getConstraintsCheckVals(constraintsStruct_scalar, ModelDCC_scalar);
end

% running the estimation of the scalar DCC model and recording into the
% resulting ModelOutLHD structure
[ModelOutLHD.outStruct_scalar, ModelOutLHD.effortStat_scalar] = ...
    estimateDCC(Sample, ModelDCC_scalar, ModelDCC_scalar, ifShowPlots);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DCC rank deficient r = 1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% estimation of the rank 1 deficient DCC model

% setting the rank value
r = 1;
% dimension of the parameter subspace for the rank deficient model
Nstar = n*r - 0.5 * r * (r - 1);

% using the estimated scalar DCC parameters as the initial ones for the
% rank one deficient models
% the user is free to change this initialization provided that the initial
% values satisfy the parameter constraints for the rank deficient DCC model
ModelDCC_rank_def1 = struct('a', ones(Nstar, 1) * sqrt(ModelOutLHD.outStruct_scalar.a), ...
    'b', ones(Nstar, 1) * sqrt(ModelOutLHD.outStruct_scalar.b),  ...
    'S', S0, 'Q', S0, 'modelType', 'rankDeficient', 'n', n, 'N', N, 'rank', r, 'K', K, 'Nstar', Nstar);

% check whether the parameter constraints are satisfied for the chosen
% initial values
constraintsStruct_rank_def = getConstraintVars(ModelDCC_rank_def1);

arrayOfConstraintsVals_rank_def = getConstraintsCheckVals(constraintsStruct_rank_def, ModelDCC_rank_def1);

% replace the initial values by the random ones if the constraints are not
% satisfied
% the initial random values are generated until the ones that satisfy the
% constraints are found
while (min(arrayOfConstraintsVals_rank_def) <= epsProjection_S)
    
    ModelDCC_rank_def1.a = randn(Nstar, 1);
    ModelDCC_rank_def1.b = randn(Nstar, 1);
    
    constraintsStruct_rank_def = getConstraintVars(ModelDCC_rank_def1);
    arrayOfConstraintsVals_rank_def = getConstraintsCheckVals(constraintsStruct_rank_def, ModelDCC_rank_def1);
end

% running the estimation of the rank one deficient DCC model and recording into the
% resulting ModelOutLHD structure
[ModelOutLHD.outStruct_rankdef1, ModelOutLHD.effortStat_rankdef1] = estimateDCC(Sample, ModelDCC_rank_def1, ModelDCC_rank_def1, ...
    ifShowPlots);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DCC rank deficient r = 2 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% estimation of the rank 2 deficient DCC model

% setting the rank value
r = 2;

% dimension of the parameter subspace for the rank 2 deficient model
Nstar = n*r - 0.5 * r * (r - 1);

% using the estimated scalar DCC parameters as the initial ones for the
% rank 2 deficient models
% the user is free to change this initialization provided that the initial
% values satisfy the parameter constraints for the rank deficient DCC model
ModelDCC_rank_def2 = struct('a', [sqrt(ModelOutLHD.outStruct_scalar.a) * ones(n, 1); ones(Nstar - n, 1) * 0.06], ...
    'b', [sqrt(ModelOutLHD.outStruct_scalar.b) * ones(n, 1); ones(Nstar - n, 1) * 0.05],  ...
    'S', S0, 'Q', S0, 'modelType', 'rankDeficient', 'n', n, 'N', N, 'rank', r, 'K', K, 'Nstar', Nstar);

% check whether the parameter constraints are satisfied for the chosen
% initial values
constraintsStruct_rank_def = getConstraintVars(ModelDCC_rank_def2);

arrayOfConstraintsVals_rank_def = getConstraintsCheckVals(constraintsStruct_rank_def, ModelDCC_rank_def2);

% replace the initial values by the random ones if the constraints are not
% satisfied
% the initial random values are generated until the ones that satisfy the
% constraints are found
while (min(arrayOfConstraintsVals_rank_def) <= epsProjection_S)
    
    ModelDCC_rank_def2.a = randn(Nstar, 1);
    ModelDCC_rank_def2.b = randn(Nstar, 1);
    
    constraintsStruct_rank_def = getConstraintVars(ModelDCC_rank_def2);
    arrayOfConstraintsVals_rank_def = getConstraintsCheckVals(constraintsStruct_rank_def, ModelDCC_rank_def2);
end

% running the estimation of the rank 2 deficient DCC model and recording into the
% resulting ModelOutLHD structure
[ModelOutLHD.outStruct_rankdef2, ModelOutLHD.effortStat_rankdef2] = estimateDCC(Sample, ModelDCC_rank_def2, ModelDCC_rank_def2, ...
    ifShowPlots);


%%%%%%%%%%%%%%%%%%%
% DCC hadamard %
%%%%%%%%%%%%%%%%%%%

% estimation of the Hadamard DCC model

% using the estimated scalar DCC parameters as the initial ones for the
% Hadamard DCC model
% the user is free to change this initialization provided that the initial
% values satisfy the parameter constraints for the Hadamard DCC model
hadamard_a = matr(ones(n, 1) * sqrt(ModelOutLHD.outStruct_scalar.a), n, 1);
hadamard_b = matr(ones(n, 1) * sqrt(ModelOutLHD.outStruct_scalar.b), n, 1);
% projecting the parameters A and B into the cone of the positive definite
% symmetric matrices
hadamard_a_proj = Projection_S(hadamard_a * hadamard_a', 'plus', 1e-4);
hadamard_b_proj = Projection_S(hadamard_b * hadamard_b', 'plus', 1e-4);

ModelDCC_hadamard = struct('a', vech(hadamard_a_proj), ...
    'b', vech(hadamard_b_proj), 'S', S0, ...
    'Q', S0, 'modelType', 'Hadamard', 'n', n, 'N', N);

% check whether the parameter constraints are satisfied for the chosen
% initial values
constraintsStruct_hadamard = getConstraintVars(ModelDCC_hadamard);

arrayOfConstraintsVals_hadamard = getConstraintsCheckVals(constraintsStruct_hadamard, ModelDCC_hadamard);

% replace the initial values by the random ones if the constraints are not
% satisfied
% the initial random values are generated until the ones that satisfy the
% constraints are found
while (min(arrayOfConstraintsVals_hadamard) < 5*epsProjection_S)
    
    ModelDCC_hadamard.a = 0.01 * randn(N, 1) + 0.01 * rand(N, 1);
    ModelDCC_hadamard.b = 0.01 * randn(N, 1) - 0.1 * rand(N, 1);
    matha = math(ModelDCC_hadamard.a);
    mathb = math(ModelDCC_hadamard.b);
    ModelDCC_hadamard.a = vech(matha * matha');
    ModelDCC_hadamard.b = vech(mathb * mathb');
    
    constraintsStruct_hadamard = getConstraintVars(ModelDCC_hadamard);
    arrayOfConstraintsVals_hadamard = getConstraintsCheckVals(constraintsStruct_hadamard, ModelDCC_hadamard);
end

% running the estimation of the Hadamard DCC model and recording into the
% resulting ModelOutLHD structure
[ModelOutLHD.outStruct_hadamard, ModelOutLHD.effortStat_hadamard] = estimateDCC(Sample, ModelDCC_hadamard, ModelDCC_hadamard, ...
    ifShowPlots);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Almon DCC              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% estimation of the Almon DCC model

% using the estimated scalar DCC parameters as the initial ones for the
% Almon DCC model
% the user is free to change this initialization provided that the initial
% values satisfy the parameter constraints for the Almon DCC model
ModelDCC_almon = struct('a', [sqrt(ModelOutLHD.outStruct_scalar.a) - 1, 0, 0]', ...
    'b', [sqrt(ModelOutLHD.outStruct_scalar.b) - 1, 0, 0]', ...
    'S', S0, 'Q', S0, 'modelType', 'almon', 'n', n, 'N', N);

% check whether the parameter constraints are satisfied for the chosen
% initial values
constraintsStruct_almon = getConstraintVars(ModelDCC_almon);

arrayOfConstraintsVals_almon = getConstraintsCheckVals(constraintsStruct_almon, ModelDCC_almon);

% replace the initial values by the random ones if the constraints are not
% satisfied
% the initial random values are generated until the ones that satisfy the
% constraints are found
while (min(arrayOfConstraintsVals_almon) <= epsProjection_S)
    
    ModelDCC_almon.a = [randn(1, 1), 0, 0]';
    ModelDCC_almon.b = [randn(1, 1), 0, 0]';
    
    constraintsStruct_almon = getConstraintVars(ModelDCC_almon);
    arrayOfConstraintsVals_almon = getConstraintsCheckVals(constraintsStruct_almon, ModelDCC_almon);
end

% running the estimation of the Almon DCC model and recording into the
% resulting ModelOutLHD structure
[ModelOutLHD.outStruct_almon, ModelOutLHD.effortStat_almon] = estimateDCC(Sample, ModelDCC_almon, ModelDCC_almon, ...
    ifShowPlots);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Almon Shuffle DCC     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% estimation of the Almon shuffle DCC model

% reorganize the order of the assets based on the PCI analysis: the asset
% logreturns are organized based of the magnitude of their projection on
% the first principal component
indexes = reorderToPC(Sample.series, n);

% creating a structure for estimation with the shuffled asset order
Sample_shfl = Sample;
Sample_shfl.series = Sample.series(indexes, :);
Sample_shfl.Epsilons = Sample.Epsilons(indexes, :);
Sample_shfl.sigmas = Sample.sigmas(indexes, :);

S0 = cov(Sample_shfl.Epsilons');

% using the estimated scalar DCC parameters as the initial ones for the
% Almon shuffle DCC model
% the user is free to change this initialization provided that the initial
% values satisfy the parameter constraints for the Almon shuffle DCC model
ModelDCC_almon_shfl = struct('a', [sqrt(ModelOutLHD.outStruct_scalar.a) - 1, 0, 0]',...
    'b', [sqrt(ModelOutLHD.outStruct_scalar.b) - 1, 0, 0]', ...
    'S', S0, 'Q', S0, 'modelType', 'almon', 'n', n, 'N', N);

% check whether the parameter constraints are satisfied for the chosen
% initial values
constraintsStruct_almon_shfl = getConstraintVars(ModelDCC_almon_shfl);

arrayOfConstraintsVals_almon_shfl = getConstraintsCheckVals(constraintsStruct_almon_shfl, ModelDCC_almon_shfl);

% replace the initial values by the random ones if the constraints are not
% satisfied; the initial random values are generated until the ones that satisfy the
% constraints are found
while (min(arrayOfConstraintsVals_almon_shfl) <= epsProjection_S)
    
    ModelDCC_almon_shfl.a = [randn(1, 1), 0, 0]';
    ModelDCC_almon_shfl.b = [randn(1, 1), 0, 0]';
    
    constraintsStruct_almon_shfl = getConstraintVars(ModelDCC_almon_shfl);
    arrayOfConstraintsVals_almon_shfl = getConstraintsCheckVals(constraintsStruct_almon_shfl, ModelDCC_almon_shfl);
end

% running the estimation of the Almon shuffle DCC model and recording into the
% resulting ModelOutLHD structure
[ModelOutLHD.outStruct_almon_shfl, ModelOutLHD.effortStat_almon_shfl] = estimateDCC(Sample_shfl, ...
    ModelDCC_almon_shfl, ModelDCC_almon_shfl, ifShowPlots);


% save the results into the file
save('res.mat');



