function pndmatrix = PnDmatrix(n)
  
    N = .5 * n * (n + 1);
    pndmatrix = zeros(N,N);
    
    for i = 1:n
        comp = sig(i, i, n);
        pndmatrix(comp,comp) = 1;
    end

