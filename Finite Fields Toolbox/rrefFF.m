function Arref = rrefFF(A,modulo)
% rrefFF(A,modulo) computes the row-reduced echelon form on matrices in 
% finite fields
% 
% Graham Holt, Aug. 2026. Updated Aug. 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% rrefFF(A,modulo)
% Arref = rrefFF(A,modulo)
% 
%% Description
% rrefFF(A,modulo) performs RREF reduction on the matrix "A" in the finite
% field "modulo"

Arref = mod(A,modulo);
[n,m] = size(A);

for k = 1:min(n,m)
    nzRow = find(gcd(Arref(k:n,k),modulo)==1,1,'first');
    nzRow = nzRow + k - 1;
    
    if isempty(nzRow)

        Arref(k:end,k) = zeros(n-k+1,1);

    else
        
        Arref([k nzRow],k:m) = Arref([nzRow k],k:m);
        Arref(k,k:m) = mod(Arref(k,k:m)*invFF(Arref(k,k),modulo),modulo);

        for j = [1:(k-1) (k+1):n]
            Arref(j,k:m) = mod(Arref(j,k:m) - Arref(j,k)*Arref(k,k:m),modulo);
        end

    end
end