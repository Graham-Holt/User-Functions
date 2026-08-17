function C = cof(M)
% cof(M) computes cofactor matrices.
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% cof(M)
% C = cof(M)
% 
%% Description
% cof(M) returns the cofactor matrix as computed by cofactor determinants.

n = size(M,1);

for k = 1:n
    indCol = [1:(k-1) (k+1):n];
    for j = 1:n
        indRow = [1:(j-1) (j+1):n];
        sign = 1 - 2*mod(j+k,2);

        C(j,k) = sign*det(M(indRow,indCol));
    end
end