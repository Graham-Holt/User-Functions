function C = cof(M)

n = size(M,1);

for k = 1:n
    indCol = [1:(k-1) (k+1):n];
    for j = 1:n
        indRow = [1:(j-1) (j+1):n];
        sign = 1 - 2*mod(j+k,2);

        C(j,k) = sign*det(M(indRow,indCol));
    end
end