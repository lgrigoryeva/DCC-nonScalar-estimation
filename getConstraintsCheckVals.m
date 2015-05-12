function arrayOfVals = getConstraintsCheckVals(constraintsStruct, Model)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

% depending on the model specification the values of the associated
% constraints are computed for the considered parameter values
    if strcmp(Model.modelType, 'scalar')
        arrayOfVals = zeros(3, 1);
        arrayOfVals(1) = min(constraintsStruct.F);
        arrayOfVals(2) = Model.a;
        arrayOfVals(3) = Model.b;
    elseif strcmp(Model.modelType, 'Hadamard')
        arrayOfVals = zeros(3, 1);
        arrayOfVals(1) = min(real(eig(math(Model.a))));
        arrayOfVals(2) = min(real(eig(math(Model.b))));
        arrayOfVals(3) = min(real(eig(constraintsStruct.D)));
    elseif strcmp(Model.modelType, 'rankDeficient')
        arrayOfVals = zeros(4, 1);
        arrayOfVals(1) = min(real(eig(constraintsStruct.D)));
        arrayOfVals(2) = min(constraintsStruct.I);
        arrayOfVals(3) = min(constraintsStruct.Kp);
        arrayOfVals(4) = min(constraintsStruct.Kn);
    elseif strcmp(Model.modelType, 'almon')
        arrayOfVals = zeros(2, 1);
        arrayOfVals(1) = min(real(eig(constraintsStruct.D)));
        arrayOfVals(2) = min(constraintsStruct.I);
    else
        display('wrong model type');
        return;
    end

end

