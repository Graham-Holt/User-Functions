function K = placeX(A,B,eig)
% placeX(A,B,eig) places poles for state feedback design.
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% placeY(A,B,eig)
% K = placeY(___)
% 
%% Description
% placeY(A,B,eig) numerically solves for a gain matrix to place closed-loop 
% poles in the provided linearized system.

[n,m] = size(B);
if length(eig) < n
    error('Too few eigenvalues to place.');
elseif length(eig) > n
    error('Too many eigenvalues to place');
end

K = sym('k%d%d',[m n]);

S = solve(charpoly(A - B*K) == poly(eig), symvar(K));

for k = 1:m
    for j = 1:n
        K(k,j) = S.("k"+num2str(k)+num2str(j));
    end
end
K = double(K);