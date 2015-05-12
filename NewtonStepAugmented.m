function out_h = NewtonStepAugmented(Model,...
        ModelTilde, ga, gb, Hess, ...
        L, tolGsmres, numberInnerIter)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva
    
    % run Newton-Raphson
    Gvector = getGradientLocalModel(Model,...
        ModelTilde, ga, gb, Hess, L);

    constraintsStruct = getConstraintVars(Model);
    constraintsStructTilde = getConstraintVars(ModelTilde);

    [out_h, ~] = gmres(@(h) TangentAugmented(h, Hess, Model, constraintsStruct, constraintsStructTilde, L), -Gvector, numberInnerIter, tolGsmres);






