function z = sig(i, j, n)
    %Given a position (i,j) in a matrix of size n gives the corresponding
    %coordinate in its vech representation
    if i >= j
        z = (n - .5 * j) * (j - 1) + i;
    else
        z = (n - .5 * i) * (i - 1) + j;
    end
    

