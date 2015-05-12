function Gvector = getGradientLocalModel(Model,...
        ModelTilde, ga, gb, Hess, L)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % computes the gradient of the local model in a and b diredtions
    % Model - current iteration; ModelTilde - previous iteration
    Ga = ga;
    Gb = gb;

    ImHessian = Hess * [Model.a - ModelTilde.a; Model.b - ModelTilde.b];

    constraintsStructTilde = getConstraintVars(ModelTilde);
    constraintsStruct = getConstraintVars(Model);

    
    if strcmp(Model.modelType, 'scalar')
        Fin = 1./constraintsStruct.F - 1./constraintsStructTilde.F;
        sumF = Fin;
        sumA = 1/ModelTilde.a - 1/Model.a;
        sumB = 1/ModelTilde.b - 1/Model.b;
        Ga = Ga + ImHessian(1) + L * (sumA + sumF);
        Gb = Gb + ImHessian(2) + L * (sumB + sumF);
    elseif strcmp(Model.modelType, 'Hadamard')
        N = Model.N;
        ImHessiana = ImHessian(1:N);
        ImHessianb = ImHessian(N + 1:2 * N);
        sumA = mathstar(inv(math(ModelTilde.a)) - inv(math(Model.a)));
        sumB = mathstar(inv(math(ModelTilde.b)) - inv(math(Model.b)));
        sumD = mathstar((inv(constraintsStructTilde.D) - inv(constraintsStruct.D)).*Model.S);
        Ga = Ga + ImHessiana + L * (sumA - sumD);
        Gb = Gb + ImHessianb + L * (sumB - sumD);
    elseif strcmp(Model.modelType, 'rankDeficient')
        N = Model.N;
        n = Model.n;
        r = Model.rank;
        Nstar = Model.Nstar;
        ImHessiana = ImHessian(1:Nstar);
        ImHessianb = ImHessian(Nstar + 1:2 * Nstar);
        ahat = matr(Model.a, n, r);
        bhat = matr(Model.b, n, r);
        
        temp = ((inv(constraintsStructTilde.D) - inv(constraintsStruct.D)).*Model.S);
        sumDA = vecr(temp * ahat);
        sumDB = vecr(temp * bhat);
        
        Kinp = 1./constraintsStruct.Kp - 1./constraintsStructTilde.Kp;
        Kinn = - 1./constraintsStruct.Kn + 1./constraintsStructTilde.Kn;
        ei = zeros(Nstar, 1);
        sumKpA = 0;
        sumKpB = 0;
        sumKnA = 0;
        sumKnB = 0;
        for i = 1 : Nstar
            ei(i) = 1;
            sumKpA = sumKpA + Kinp(i) * ei;
            sumKpB = sumKpB + Kinp(i + Nstar) * ei;
            sumKnA = sumKnA + Kinn(i) * ei;
            sumKnB = sumKnB + Kinn(i + Nstar) * ei;
            ei(i) = 0;
        end
        
        sumI = 1./constraintsStructTilde.I - 1./constraintsStruct.I;
        
        resI = zeros(Nstar, 1);
        resI(1) = sumI(1);
        Ga = Ga + ImHessiana + L * (2 * (- sumDA) + resI + sumKpA + sumKnA);
        resI(1) = sumI(2);
        Gb = Gb + ImHessianb + L * (2 * (- sumDB) + resI + sumKpB + sumKnB);
    elseif strcmp(Model.modelType, 'almon')
        ImHessiana = ImHessian(1:3);
        ImHessianb = ImHessian(4:6);
        n = Model.n;
        ahat = almonFunction(n, Model.a);
        bhat = almonFunction(n, Model.b);
        la = ahat - Model.a(1);
        lb = bhat - Model.b(1);
        ind1 = (1:n)';
        ind2 = ind1.^2;
        i_n = ones(n, 1);
        DaT = [i_n ind1.*la ind2.*la]';
        DbT = [i_n ind1.*lb ind2.*lb]';
        N = Model.N;
        temp = ((inv(constraintsStructTilde.D) - inv(constraintsStruct.D)).*Model.S);
        sumDA = DaT * (temp * ahat);
        sumDB = DbT * (temp * bhat);        
        
        sumI = 1./constraintsStructTilde.I - 1./constraintsStruct.I;
        e1 = zeros(n, 1);
        e1(1) = 1;
        DaTe1 = DaT * e1;
        DbTe1 = DbT * e1;
        resI = sumI(1) * DaTe1;
        Ga = Ga + ImHessiana + L * (2 * ( - sumDA) + resI);
        resI = sumI(2) * DbTe1;
        Gb = Gb + ImHessianb + L * (2 * ( - sumDB) + resI);
    end

    
    Gvector = [Ga; Gb];

end

