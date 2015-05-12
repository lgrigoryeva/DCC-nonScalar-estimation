function l = elimination(n)
% Author: Juan-Pablo Ortega

    % Returns  elimination matrix of size n

    l = zeros(n * (n + 1) / 2, n^2);
    a = zeros(n, n);
    for i = 1:n
        for j = 1:n
            a(i, j) = 1;
            l(:,(j - 1) * n + i) = vech(a);
            a(i, j) = 0;
        end
    end