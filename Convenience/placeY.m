function K = placeY(A,B,C,D,eig)
% placeY(A,B,C,D,eig) places poles for output feedback design.
% 
% Graham Holt, August 2026. Updated September 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% placeY(A,B,C,D,eig)
% K = placeY(___)
% 
%% Description
% placeY(A,B,C,D,eig) numerically solves for a gain matrix to place closed-
% loop poles in the provided linearized system.

[n,m] = size(B);
[p,~] = size(C);
if length(eig) < 2*n-m-p+1
    error('Too few eigenvalues to place.');
end
s = n-m-p+1;
eig = eig(1:(n+s));

A = blkdiag(A,zeros(s));
B = blkdiag(B,eye(s));
C = blkdiag(C,eye(s));
D = blkdiag(D,zeros(s));

K = sym('k%d%d',[m p]+s);

S = solve(charpoly(A - B*K/(eye(p+s) + D*K)*C) == poly(eig), symvar(K));

for k = 1:(m+s)
    for j = 1:(p+s)
        K(k,j,:) = reshape(S.("k"+num2str(k)+num2str(j)),1,1,[]);
    end
end
K = double(K);
[~,I] = mink(reshape(pagenorm(K),1,[]),1);
K = K(:,:,I);