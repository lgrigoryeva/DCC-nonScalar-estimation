function S3 = StructPlusStruct(S1, S2, fieldNamesList)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva
    
    % provides the sum of two structures when the fields are provided in 
    % fieldNamesList
    S3 = S1;
    for i = 1:length(fieldNamesList)
        fieldName = fieldNamesList(i);
        if isfield(S1, char(fieldName)) && isfield(S2, char(fieldName))
            S3 = setfield(S3, char(fieldName), S1.(char(fieldName)) + S2.(char(fieldName)));
        else
            display('problem');
        end
    end

end

