function [output_vector] = almonFunction(n, theta)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % this function constructs the n Almon lag function values for the
    % parameters theta

    % theta is the column or row vector of the length 3, where
    % theta(1) is the intercept
    output_vector = zeros(n, 1);
    theta_vector = theta(:);
    
    for i = 1:n
        output_vector(i) = exp(theta_vector(2) * i + theta_vector(3) * i^2) + theta_vector(1);
    end
    
    
end

