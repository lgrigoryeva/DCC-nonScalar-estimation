function [ModelOut, effortStatStruct] = estimateDCC(Sample, ModelReal, ModelInitial, ...
                                                    ifShowPlots, toleranceParamsArray)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % Likelihood optimization
    % time the start of the estimation 
    tic;

    if ifShowPlots
        clf;
    end

    n = ModelReal.n;
    N = ModelReal.N;
    
    if strcmp(ModelReal.modelType,'Hadamard')
        dimParam = N;
    elseif strcmp(ModelReal.modelType,'scalar')
        dimParam = 1;
    elseif strcmp(ModelReal.modelType,'rankDeficient')
        dimParam = ModelReal.Nstar;
    elseif strcmp(ModelReal.modelType,'almon')
        dimParam = 3;
    end
    counter = 0;
    % if the toleranceParamsArray with the options for the estimation is
    % not provided, then we fill in the default values
    if nargin < 5
        TolFun = 1e-6;
        TolX = 1e-5;
        TolNewton = 1e-6;
        numberInnerIter = dimParam * 2;
        NewtonLoopsNum = 100;
        tolGsmres = 1e-10;
        epsProjection_S = 0;
        maxTotalLoops = 2000;
    % otherwise we use the parameters supplied in toleranceParamsArray
    else
        TolFun = toleranceParamsArray(1);
        TolX = toleranceParamsArray(2);
        TolNewton = toleranceParamsArray(3);
        numberInnerIter = toleranceParamsArray(4);
        NewtonLoopsNum = toleranceParamsArray(5);
        tolGsmres = toleranceParamsArray(6);
        epsProjection_S = toleranceParamsArray(7);
        maxTotalLoops = toleranceParamsArray(8);
    end

    % initial values of the loglikelihood and its gradient for the given
    % model specification and the provided sample
    [fVal, ga, gb] = dccLikelihoodAndScore(ModelInitial, Sample);

    % initial values of the Bregman penalization strengths
    L = norm([ga; gb]);
    
    if strcmp(ModelReal.modelType, 'Hadamard')
        % we reduce the value of L for Hadamard: a matter of choice and
        % tuning
        L = L/10;
        TolX = 1e-6;
    elseif strcmp(ModelReal.modelType, 'rankDeficient')
        % we reduce the value of L for rankDeficient: a matter of choice and
        % tuning
        L = L/10;
    elseif strcmp(ModelReal.modelType, 'almon')
        % we increase the value of L for Almon: a matter of choice and
        % tuning
        L = L * 1000;
        % we change the values of TolX and TolNewton for Almon: a matter of 
        % choice and tuning
        TolX = 1e-6;
        TolNewton = 1e-9;
    end
    
    % in case the plotting is switched on
    if ifShowPlots
        fValGoodLoops = [];
        ErrorA = [];
        ErrorB = [];
        R = [];
    end

    % Starting point of the BFGS Hessian
    Hess = L * eye(2 * dimParam);
    goodLoops = 1;
    totalLoops = 1;
    TolDiffFun = [];
    TolDiffX = [];
    TolDiffFun(1) = 2 * TolFun;
    TolDiffX(1) = 2 * TolX;
    ModelCur = ModelInitial;

    ModelH = ModelInitial;

    % run the estimation until the corresponding tolerances are achieved
    % and the maximum of the iterations number is not reached
    while (TolDiffFun(goodLoops) > TolFun || TolDiffX(goodLoops) > TolX) ...
            && totalLoops < maxTotalLoops
        if ifShowPlots
            fValGoodLoops(goodLoops) = fVal;
        end

        ftilde = fVal;
        gatilde = ga;
        gbtilde = gb;
        ModelTilde = ModelCur;
        newtonLoops = 0;

        flagEnd = 0;
        flagBadNewton = 0;
        flagBadNewtonConstr = 0;
        while ~flagEnd
            newtonLoops = newtonLoops + 1;

            h = NewtonStepAugmented(ModelCur, ...
                ModelTilde, ga, gb, Hess, L, tolGsmres, numberInnerIter);
            if isnan(h)
                flagEnd = 1;
                flagBadNewtonConstr = 1;
                break;
            end
            ModelH.a = h(1:dimParam);
            ModelH.b = h(dimParam + 1:2 * dimParam);


            ModelParPlusHpar = StructPlusStruct(ModelH, ModelCur, {'a', 'b'});
            constraintsStructParPlusHpar = getConstraintVars(ModelParPlusHpar);
            arrayOfConstraintsVals = getConstraintsCheckVals(constraintsStructParPlusHpar, ModelParPlusHpar);
            % Division by 2
            numIter = 0;

            while (min(arrayOfConstraintsVals) <= epsProjection_S)
                numIter = numIter + 1;

                ModelH.a = ModelH.a / 2;
                ModelH.b = ModelH.b / 2;
                ModelParPlusHpar = struct();
                ModelParPlusHpar = StructPlusStruct(ModelH, ModelCur, {'a', 'b'});
                constraintsStructParPlusHpar = getConstraintVars(ModelParPlusHpar);
                arrayOfConstraintsVals = getConstraintsCheckVals(constraintsStructParPlusHpar, ModelParPlusHpar);

            end


            ModelCur = StructPlusStruct(ModelH, ModelCur, {'a', 'b'});


            Gvector = getGradientLocalModel(ModelCur, ModelTilde, ga, gb, Hess, L);

            tolNewton = norm(Gvector);

            if newtonLoops > NewtonLoopsNum
                constraintsStruct = getConstraintVars(ModelCur);
                arrayOfConstraintsVals = getConstraintsCheckVals(constraintsStruct, ModelCur);
                if min(arrayOfConstraintsVals) <= epsProjection_S
                    flagBadNewtonConstr = 1;
                    counter = counter + 1;
                else
                    flagBadNewtonConstr = 0;
                end
                flagEnd = 1;
                flagBadNewton = 1;
            end
            if tolNewton < TolNewton
                flagEnd = 1;
            end

        end
        display(newtonLoops);



        % Trust region evaluation
        sA = ModelCur.a - ModelTilde.a;
        sB = ModelCur.b - ModelTilde.b;

        % Trust region: a matter of tuning
        if strcmp(ModelReal.modelType, 'rankDeficient')
            expectedDecrease = gatilde' * sA + gbtilde' * sB ;
            condExpand = 'rho > 0.9 &&   rho < 2';% && ~flagBadNewton';
            condShrink = 'rho < 0.01 || rho > 2';
        elseif strcmp(ModelReal.modelType, 'almon')
            expectedDecrease = gatilde' * sA + gbtilde' * sB ;
            condExpand = 'rho >= 0.8  && rho < 1.5';% && ~flagBadNewton';
            %condExpand = 'rho > 0.9  && rho < 2 && ~flagBadNewton';
            condShrink = 'rho < 0.01 || rho > 1.5';
            %condShrink = 'rho < 0.01 || rho > 2';
        elseif strcmp(ModelReal.modelType, 'Hadamard')
            expectedDecrease = gatilde' * sA + gbtilde' * sB;
            condExpand = 'rho > 0.8 && rho < 2';% && ~flagBadNewton';%0.9
            condShrink = 'rho < 0.01 || rho > 2';%no rho>2
        else
            expectedDecrease = gatilde' * sA + gbtilde' * sB;
            condExpand = 'rho >= 0.8';% && ~flagBadNewton';
            condShrink = 'rho < 0.01';
        end


        fVal = dccLikelihood(ModelCur, Sample);

        rho = (fVal - ftilde) / expectedDecrease;

        if ifShowPlots
            if ~isnan(rho) && ~imag(rho)
                R(totalLoops) = rho;
                Lvals(totalLoops) = L;
            end
        end
        flag = 1;
        display(counter);
        % change the value of the penalization strength depending on the
        % outcome of the iteration
        if isnan(rho) || eval(condShrink) || imag(rho)% || flagBadNewtonConstr %|| flagBadNewton
            L = 2 * L;
            flag = 0;
        elseif eval(condExpand) && ~flagBadNewtonConstr %&& ~flagBadNewton
            L = .5 * L;
            counter = 0;
        end
        if counter > 40
            break;
        end

        if flag && ~flagBadNewtonConstr
            [fVal, ga, gb] = dccLikelihoodAndScore(ModelCur, Sample);
            Astruct.ga(goodLoops).ga = ga;
            Bstruct.gb(goodLoops).gb = gb;

            if ifShowPlots
                ErrorA(goodLoops) = norm(ModelCur.a - ModelReal.a);
                ErrorB(goodLoops) = norm(ModelCur.b - ModelReal.b);
                Astruct.ErrorA(goodLoops) = ErrorA(goodLoops);
                Astruct.Aval(goodLoops).a = ModelCur.a;
                Bstruct.ErrorB(goodLoops) = ErrorB(goodLoops);
                Bstruct.Bval(goodLoops).b = ModelCur.b;

            end
            % update the Hessian
            s = [sA; sB];
            yA = ga - gatilde;
            yB = gb - gbtilde;
            y = [yA; yB];
            yT = y';
            Hesstimess = Hess * s;
            sT = s';
            Hess = Hess + (y * yT)/(yT * s) - (Hesstimess * Hesstimess')/(sT * Hesstimess);
            goodLoops = goodLoops + 1;
            TolDiffFun(goodLoops) = abs((ftilde - fVal) / ftilde);
            TolDiffX(goodLoops) = norm(ModelCur.a - ModelTilde.a)/norm(ModelTilde.a) + ...
                norm(ModelCur.b - ModelTilde.b)/norm(ModelTilde.b);
            display(TolDiffFun(goodLoops));
            display(TolDiffX(goodLoops));
        else
            fVal = ftilde;
            ga = gatilde;
            gb = gbtilde;
            ModelCur = ModelTilde;
            counter = counter + 1;
        end
        totalLoops = totalLoops + 1;
        if ifShowPlots
            figure(1)
            subplot(4,1,1)
            plot(ErrorA)
            subplot(4,1,2)
            plot(ErrorB)
            subplot(4, 1, 3)
            plot(R)
            subplot(4, 1, 4)
            plot(fValGoodLoops)
            drawnow
        end


    end

    ModelOut = ModelCur;

    elapsedTime = toc;
    % save all the details into the structure
    if ~ifShowPlots
        fValGoodLoops = [];
        R = [];
        Lvals = [];
    end
    effortStatStruct = struct('A', Astruct, 'B', Bstruct, 'fValGoodLoops', fValGoodLoops, 'TolX', TolDiffX, 'TolFun', TolDiffFun, ...
        'R', R, 'L', Lvals, 'elapsedTime', elapsedTime, 'totalLoops', totalLoops, 'goodLoops', goodLoops, 'ifTotalLoopsExceeded', 0);
    if totalLoops >= maxTotalLoops
        effortStatStruct.ifTotalLoopsExceeded = 1;
    end
