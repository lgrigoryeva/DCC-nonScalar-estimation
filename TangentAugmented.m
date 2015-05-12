function Tvector = TangentAugmented(h, Hess, Model, constraintsStruct, constraintsStructTilde, L)       
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva
    
    % compute the Jacobian of the gradient of the local model at the
    % current iteration
    if strcmp(Model.modelType, 'scalar')
        Fin2 = 1./(constraintsStruct.F.^2);
        Ha = h(1);
        Hb = h(2);
        ImHessian = Hess * [Ha; Hb];
        sumA = Ha/Model.a^2;
        sumB = Hb/Model.b^2;
        sumF = Fin2 * (Ha + Hb);
        Ta = ImHessian(1) + L * (sumA + sumF);
        Tb = ImHessian(2) + L * (sumB + sumF);
    elseif strcmp(Model.modelType, 'Hadamard')
        N = Model.N;
        Ha = h(1:N);
        Hb = h(N + 1:2 * N);

        ImHessian = Hess * [Ha; Hb];

        ImHessiana = ImHessian(1:N);
        ImHessianb = ImHessian(N + 1:2 * N);

        
        invA = inv(math(Model.a));
        invB = inv(math(Model.b));
        invD = inv(constraintsStruct.D);
        
        sumA = mathstar(invA * math(Ha) * invA);
        sumB = mathstar(invB * math(Hb) * invB);
        sumD = mathstar((invD * (math(Ha + Hb).*Model.S) * invD).*Model.S);
        
        Ta = ImHessiana + L * (sumA + sumD);
        Tb = ImHessianb + L * (sumB + sumD);
    elseif strcmp(Model.modelType, 'rankDeficient')

        N = Model.N;
        n = Model.n;
        r = Model.rank;
        Nstar = Model.Nstar;
        Ha = h(1:Nstar);
        Hb = h(Nstar + 1:2 * Nstar);

        ImHessian = Hess * [Ha; Hb];

        ImHessiana = ImHessian(1:Nstar);
        ImHessianb = ImHessian(Nstar + 1:2 * Nstar);
                
        ahat = matr(Model.a, n, r);
        bhat = matr(Model.b, n, r);
        Hahat = matr(Ha, n, r);
        Hbhat = matr(Hb, n, r);
        
        Htheta = Hahat * ahat' + ahat * Hahat' + Hbhat * bhat' + bhat * Hbhat';

        invD = inv(constraintsStruct.D);
        invDtilde = inv(constraintsStructTilde.D);
        temp1 = (invD * (Htheta.*Model.S) * invD).*Model.S;
        temp2 = (invDtilde - invD).*Model.S;
        
        sumDA = vecr(temp1 * ahat - temp2 * Hahat);
        sumDB = vecr(temp1 * bhat - temp2 * Hbhat);
        
        Kin2p = [Ha; Hb]./(constraintsStruct.Kp.^2);
        Kin2n = [Ha; Hb]./(constraintsStruct.Kn.^2);
        sumKpA = 0;
        sumKpB = 0;
        sumKnA = 0;
        sumKnB = 0;
        ei = zeros(Nstar, 1);
        for i = 1:Nstar
            ei(i) = 1;
            sumKpA = sumKpA + Kin2p(i) * ei;
            sumKpB = sumKpB + Kin2p(i + n) * ei;
            sumKnA = sumKnA + Kin2n(i) * ei;
            sumKnB = sumKnB + Kin2n(i + n) * ei;
            ei(i) = 0;
        end
        
        sumI = [Ha(1); Hb(1)]./constraintsStruct.I.^2;
        resI = zeros(Nstar, 1);
        resI(1) = sumI(1);
        Ta = ImHessiana + L * (2 * (sumDA) + resI + sumKpA + sumKnA);
        resI(1) = sumI(2);
        Tb = ImHessianb + L * (2 * (sumDB) + resI + sumKpB + sumKnB);
    elseif strcmp(Model.modelType, 'almon')

        n = Model.n;

        Ha = h(1:3);
        Hb = h(4:6);

        ImHessian = Hess * [Ha; Hb];

        ImHessiana = ImHessian(1:3);
        ImHessianb = ImHessian(4:6);
        
        
        ahat = almonFunction(n, Model.a);
        bhat = almonFunction(n, Model.b);
        la = ahat - Model.a(1);
        lb = bhat - Model.b(1);
        ind1 = (1:n)';
        ind2 = ind1.^2;
        i_n = ones(n, 1);
        Da = [i_n ind1.*la ind2.*la];
        Db = [i_n ind1.*lb ind2.*lb];
        DaT = Da';
        DbT = Db';
        
        Hahat = Da * Ha;
        Hbhat = Db * Hb;
        
        Htheta = Hahat * ahat' + ahat * Hahat' + Hbhat * bhat' + bhat * Hbhat';
        %%%%%
        zero_n = zeros(n, 1);
        Da1 = [zero_n ind1.*la ind2.*la];
        Db1 = [zero_n ind1.*lb ind2.*lb];
        Hahat1 = Da1 * Ha;
        Hbhat1 = Db1 * Hb;
        
        D2a = [zero_n ind1.*Hahat1 ind2.*Hahat1];
        D2b = [zero_n ind1.*Hbhat1 ind2.*Hbhat1];
        
        H2ahatT = (D2a * Ha)';
        H2bhatT = (D2b * Hb)';

        
        invD = inv(constraintsStruct.D);
        invDtilde = inv(constraintsStructTilde.D);
        temp1 = - (invD * (Htheta.*Model.S) * invD).*Model.S;
        
        temp2 = (invDtilde - invD).*Model.S;
        
        sumDA = H2ahatT * temp2 * ahat + DaT * temp1 * ahat + DaT * temp2 * Hahat;
        sumDB = H2bhatT * temp2 * bhat + DbT * temp1 * bhat + DbT * temp2 * Hbhat;
        e1 = zeros(n, 1);
        e1(1) = 1;
        DaTe1 = DaT * e1;
        DbTe1 = DbT * e1;
        D2aTe1 = H2ahatT * e1;
        D2bTe1 = H2bhatT * e1;
        sumI1 = [Hahat(1); Hbhat(1)]./constraintsStruct.I.^2;
        sumI2 = (1./constraintsStructTilde.I - 1./constraintsStruct.I);
        
        resI = sumI1(1) * DaTe1 + sumI2(1) * D2aTe1;
        Ta = ImHessiana + L * (2 * (- sumDA) + resI);
        resI = sumI1(2) * DbTe1 + sumI2(2) * D2bTe1;
        Tb = ImHessianb + L * (2 * (- sumDB) + resI);
            
    end

    Tvector = [Ta; Tb];
    
end 

