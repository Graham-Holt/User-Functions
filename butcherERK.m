function S = butcherERK(p)
% butcherERK(p) generates the order conditions for a stable, p-order,
% explicit Runge-Kutta Method.
% 
% Graham Holt, June 2026. Updated June 2026.
% Embry-Riddle Aeronautical University
% 
%% Syntax
% butcherERK(p)
% S = butcherERK(___)
% 
%% Description
% butcherERK(p) generates the order conditions for a stable, p-order,
% explicit Runge-Kutta Method.

tau = 0;
for k = 1:p
    tau = tau + numRootTree(k);
end

s = ceil(0.5*(sqrt(1 + 8*tau) - 1)); Phi = {};
for k = 1:p
    Phi{k} = cell2mat(rootTree(k,s));
end
Phi = cell2mat(Phi); gamma = Phi(end,:).'; Phi = Phi(1:(end-1),:).';

A = sym('a%d%d',[s s]); b = sym('b%d',[s 1]);  c = sym('c%d',[s 1]);

S = [Phi*b==gamma; c==tril(A,-1)*ones(s,1)];
disp(['Order conditions require ',num2str(s),' steps and ',num2str(0.5*s*(s+1)-size(Phi,1)),' free parameters']);

end

function T = rootTree(p,s)
% rootTree(p,s) generates the elementary weights and reciprocal density of
% a p-order rooted tree summed over "s" steps.
% 
% Graham Holt, June 2026. Updated June 2026.
% Embry-Riddle Aeronautical University
% 
%% Syntax
% rootTree(p,s)
% S = rootTree(___)
% 
%% Description
% rootTree(p,s) returns the elementary weights and reciprocal density of
% a p-order rooted tree summed over "s" steps as a cell array.

if s<p
    error('Number of steps must be greater than or equal to the order');
end
if p==1
    T = {ones(s+1,1)};
    return;
end

A = sym('a%d%d',[s s]); A = blkdiag(tril(A,-1),1);
K = intpartition(p-1); k = 1;
for j = 1:length(K)
    Phi = {}; M = [];
    for i = 1:length(K{j})
        Phi{i} = rootTree(K{j}(i),s);
        M(i) = size(Phi{i},2);
    end
    I = ones(1,length(Phi));
    while all(I<=M)
        T{k} = [ones(s,1); 1/p];
        for m = 1:length(Phi)
            T{k} = T{k}.*(A*Phi{m}(:,I(m)));
        end

        k = k + 1;
        if all(I==M)
            break;
        end

        I(end) = I(end)+1;
        while any(I>M)
            I(find(I>M,1,'last')-1) = I(find(I>M,1,'last')-1)+1;
            I(find(I>M,1,'last')) = 1;
        end
    end
end

end

function T = numRootTree(p)
% numRootTree(p) counts the number of rooted trees of order "p".
% 
% Graham Holt, June 2026. Updated June 2026.
% Embry-Riddle Aeronautical University
% 
%% Syntax
% numRootTree(p)
% T = numRootTree(___)
% 
%% Description
% numRootTree(p) returns the number of rooted trees of order "p".

if p==1
    T = 1; return;
end

T = 0; P = intpartition(p-1);
for k = 1:length(P)
    tau = 1;
    for j = 1:length(P{k})
        tau = tau*numRootTree(P{k}(j));
    end
    T = T + tau;
end

end


function P = intpartition(n)
% intpartition(n) generates the integer partitions of "n".
% 
% Graham Holt, June 2026. Updated June 2026.
% Embry-Riddle Aeronautical University
% 
%% Syntax
% intpartition(n)
% P = intpartition(n)
% 
%% Description
% intpartition(n) returns a cell array with the integer partitions of "n".

p = zeros(1,n); p(1) = n;
P = {}; k = 1;

while true
    P{end+1} = p(1:k);
    if k==n
        return;
    end

    r = k - find(p>1,1,'last') + 1;
    k = find(p>1,1,'last');
    p(k) = p(k) - 1;

    while r>p(k)
        p(k+1) = p(k);
        r = r - p(k);
        k = k + 1;
    end

    p(k+1) = r;
    k = k + 1;
end

end