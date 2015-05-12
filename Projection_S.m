function A_projected = Projection_S(A, plus_or_minus, epsProjection)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

    % Projects a n x r matrix (r<=n) onto the cone of full rank matrices with 
    % positive or negative singular eigenvalues.
    % plus_or_minus can take the values 'plus' or
    % 'minus'. Notice that this projection does not respect the space of
    % n-symmetric matrices.
    [V,D] = eig(A);
    [n,r] = size(D);
    d = diag(D);
    
    if strcmpi(plus_or_minus, 'plus')
        d1 = real(max(d, epsProjection));
    else
        d1 = real(min(d, - epsProjection));
    end
    
    if d == d1
        A_projected = A;
    else
        D1= zeros(n, r);
        for i = 1:r
           D1(i,i) = d1(i); 
        end
        A_projected = V * D1 * V';
    end
    
