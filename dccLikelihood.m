function likelihood = dccLikelihood(Model, Sample)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

% this function computes MINUS the likelihood
    
    % Model structure contains:
    % modelType - DCC model specification name (for instance: 'Hadamard',
    % 'scalar', 'rankDeficient', 'almon')
    % a, b are the intrinsic parameters of the DCC model
    % Q - the initial value Q0
    % S - cov(EE')
    % n - dimension
    
    % Sample structure contains:
    % series - the time series/observations
    % H0 - the initial covariance matrix
    % z0 - the initial time series value (n-dimensional) 
    % sigmas - the associated sonditional deviations computed at the 1st
    % stage of the DCC model estimation by fitting the univariate
    % GARCH(1,1) model to the observed series
  
z = Sample.series;
H0 = Sample.H0;
z0 = Sample.z0;
epsilons0 = sqrtm(H0) \ z0; 
sigmas = Sample.sigmas;

% conputing the standardized returns
epsilons = z./sigmas;

[~, T] = size(z);
n = Model.n;
i_n = ones(n, 1);
iit = i_n * i_n';

% recovering the parameter matrices depending on the instrinsic parameters
% of the DCC model considered
if strcmp(Model.modelType, 'Hadamard')
    A = math(Model.a);
    B = math(Model.b);
elseif strcmp(Model.modelType, 'scalar')
    A = Model.a * iit;
    B = Model.b * iit;
elseif strcmp(Model.modelType, 'rankDeficient')
    r = Model.rank;
    A = matr(Model.a, n, r) * matr(Model.a, n, r)';
    B = matr(Model.b, n, r) * matr(Model.b, n, r)';
elseif strcmp(Model.modelType, 'almon')
    atilde = almonFunction(n, Model.a);
	btilde = almonFunction(n, Model.b);
	A = atilde * atilde';
	B = btilde * btilde';
else
    display('wrong model type');
    return;
end
S = Model.S;
Q = Model.Q;


% the first iteration
Q0 = Q;
Q1 = (iit - A - B).* S + A.*(epsilons0 * epsilons0') + B.*Q0;
Q1_ast_inv = diag(1./sqrt(diag(Q1)));
D1 = (diag(sigmas(:, 1)));
H1 = D1 * Q1_ast_inv * Q1 * Q1_ast_inv * D1;
% loglikelihood value
l = -.5 * ((T * n) * log(2 * pi) + log(det(H1)) + z(:, 1)' * inv(H1) * z(:, 1));
Q_t = Q1;

% the next iterations
for t = 2:T
    EE = epsilons(:, t - 1) * epsilons(:, t - 1)';    
    
    Q_t = (iit - A - B).* S + A.*EE + B.*Q_t;
    Q_t_ast_inv = diag(1./sqrt(diag(Q_t)));
    
    D_t = diag((sigmas(:, t)));
    H_t = D_t * Q_t_ast_inv * Q_t * Q_t_ast_inv * D_t;

    % loglikelihood value
    l = l - 0.5 * (log(det(H_t)) + z(:, t)' * inv(H_t) * z(:, t));   
end

% The averaged likelihood with a minus sign in order to
% convert into a minimization problem
likelihood = - l/T;



