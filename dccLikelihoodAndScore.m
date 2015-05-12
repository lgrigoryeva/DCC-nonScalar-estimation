function [likelihood, gradienta, gradientb] = dccLikelihoodAndScore(Model, Sample)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % this function computes MINUS the likelihood and the gradients in a
    % and b directions
    
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
    Ln = elimination(n);
    Dn = duplication(n);
    PnD = PnDmatrix(n);
    prodPnDLn = PnD * Ln;

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
        la = atilde - Model.a(1);
        lb = btilde - Model.b(1);
        A = atilde * atilde';
        B = btilde * btilde';
        ind1 = [1:n]';
        ind2 = ind1.^2;
        Da = [i_n ind1.*la ind2.*la];
        Db = [i_n ind1.*lb ind2.*lb];
    else
        display('wrong model type');
        return;
    end
    S = Model.S;
    Q = Model.Q;
    diagVechB = diag(vech(B));

    % the first iteration
    Q0 = Q;
    Q1 = (iit - A - B).* S + A.*(epsilons0 * epsilons0') + B.*Q0;
    Q1_ast_inv = diag(1./sqrt(diag(Q1)));
    D1 = (diag(sigmas(:, 1)));
    H1 = D1 * Q1_ast_inv * Q1 * Q1_ast_inv * D1;
    H1_inv = inv(H1);
    z1 = z(:, 1);
    l = -.5 * ((T * n) * log(2 * pi) + log(det(H1)) + z1' * H1_inv * z1);
    dHl = -.5 * (H1_inv - (H1_inv * z1) * (z1' * H1_inv));

    EE = epsilons0 * epsilons0';
    A1 = diag(vech(EE - S));
    B1 = diag(vech(Q0 - S));
    iQastD = Q1_ast_inv * D1;
    i2QastD = Q1_ast_inv * iQastD;
    temp = (Ln * kron(iQastD, iQastD) - 1/2 * prodPnDLn * (kron(Q1_ast_inv * Q1 * iQastD, i2QastD) + ...
        kron(iQastD, Q1_ast_inv * Q1_ast_inv * Q1 * iQastD))) * Dn;
    ThetaHA1 = A1 * temp;
    ThetaHB1 = B1 * temp;

    dA1 = ThetaHA1 * vech(dHl);
    dB1 = ThetaHB1 * vech(dHl);

    Q_tm1 = Q1;
    A_tm1 = A1;
    B_tm1 = B1;
    dA = dA1;
    dB = dB1;

    % the next iterations
    for t = 2:T
        EE = epsilons(:, t - 1) * epsilons(:, t - 1)';    

        Q_t = (iit - A - B).* S + A.*EE + B.*Q_tm1;
        Q_t_ast_inv = diag(1./sqrt(diag(Q_t)));

        if (min(diag(Q_t)) < 0)
            display('negative Q in likelihood');
        end
        D_t = diag((sigmas(:, t)));
        H_t = D_t * Q_t_ast_inv * Q_t * Q_t_ast_inv * D_t;


        H_t_inv = inv(H_t);

        z_t = z(:, t);
        l = l - 0.5 * (log(det(H_t)) + z_t' * H_t_inv * z_t);
        dHl = - .5 * (H_t_inv - (H_t_inv * z_t) * (z_t' * H_t_inv));


        A_t = diag(vech(EE - S)) + A_tm1 * diagVechB;
        B_t = diag(vech(Q_tm1 - S)) + B_tm1 * diagVechB;

        iQastD = Q_t_ast_inv * D_t;
        i2QastD = Q_t_ast_inv * iQastD;
        temp1 = Q_t_ast_inv * Q_t * iQastD;
        temp2 = (Ln * kron(iQastD, iQastD) - 1/2 * prodPnDLn * (kron(temp1, i2QastD) + ...
            kron(iQastD, Q_t_ast_inv * temp1))) * Dn;
        ThetaHA_t = A_t * temp2;
        ThetaHB_t = B_t * temp2;

        dA = dA + ThetaHA_t * vech(dHl);
        dB = dB + ThetaHB_t * vech(dHl);

        Q_tm1 = Q_t;
        A_tm1 = A_t;
        B_tm1 = B_t;

    end

    % The averaged likelihood and gradients with a minus sign in order to
    % convert into a minimization problem
    likelihood = - l/T;
    gradientA = - dA/T;
    gradientB = - dB/T;
    if strcmp(Model.modelType, 'Hadamard')
        gradienta = 2 * gradientA - vech(diag(diag(math(gradientA))));
        gradientb = 2 * gradientB - vech(diag(diag(math(gradientB))));
    elseif strcmp(Model.modelType, 'scalar')
        gradienta = trace(math(gradientA) * iit);
        gradientb = trace(math(gradientB) * iit);
    elseif strcmp(Model.modelType, 'rankDeficient')
        gradienta = 2 * vecr(math(gradientA) * matr(Model.a, n, r));
        gradientb = 2 * vecr(math(gradientB) * matr(Model.b, n, r));
    elseif strcmp(Model.modelType, 'almon')
        gradienta = 2 * (Da' * math(gradientA) * atilde);
        gradientb = 2 * (Db' * math(gradientB) * btilde);
    end


