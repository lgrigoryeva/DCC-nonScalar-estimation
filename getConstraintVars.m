function [constraintsStruct] = getConstraintVars(modelParamsStruct)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % modelParamsStruct contains:
    % modelType - DCC model specification name (for instance: 'Hadamard',
    % 'scalar', 'rankDeficient', 'almon', 'almon_shfl')
    % a, b - the intrinsic parameters of the DCC model
    % Q - the initial value Q0
    % S - cov(EE')
    % n - dimension
    
    % constraintsStruct contains the values of the constraints 
    if strcmp(modelParamsStruct.modelType, 'scalar')
        constraintsStruct.F = 1 - modelParamsStruct.a - modelParamsStruct.b;    
    else
        n = modelParamsStruct.n;
        in = ones(n, 1);
        if strcmp(modelParamsStruct.modelType, 'Hadamard')
            constraintsStruct.D = (in * in' - math(modelParamsStruct.a) - math(modelParamsStruct.b)).*modelParamsStruct.S;
        elseif strcmp(modelParamsStruct.modelType, 'rankDeficient')
            r = modelParamsStruct.rank;
            Nstar = modelParamsStruct.Nstar;
            iNstar = ones(Nstar, 1);
            aatPlusbbt = matr(modelParamsStruct.a, n, r) * matr(modelParamsStruct.a, n, r)' + ...
                matr(modelParamsStruct.b, n, r) * matr(modelParamsStruct.b, n, r)';
            constraintsStruct.D = (in * in' - aatPlusbbt).*modelParamsStruct.S;
            constraintsStruct.I = [modelParamsStruct.a(1); modelParamsStruct.b(1)];
            constraintsStruct.Kp = modelParamsStruct.K * [iNstar; iNstar] - [modelParamsStruct.a; modelParamsStruct.b];
            constraintsStruct.Kn = modelParamsStruct.K * [iNstar; iNstar] + [modelParamsStruct.a; modelParamsStruct.b];
       elseif strcmp(modelParamsStruct.modelType, 'almon')
            atilde = almonFunction(n, modelParamsStruct.a);
            btilde = almonFunction(n, modelParamsStruct.b);
            aatPlusbbt = atilde * atilde' + btilde * btilde';
            constraintsStruct.D = (in * in' - aatPlusbbt).*modelParamsStruct.S;
            constraintsStruct.I = [atilde(1); btilde(1)];
        else
            display('the wrong model type');
            return;
        end
    end
end

