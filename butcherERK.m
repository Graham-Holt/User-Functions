function [S,BT] = butcherERK(p,method,alpha)

s = p; Phi = {};
for k = 1:p
    Phi{k} = cell2mat(rootTree(k,s));
end
Phi = cell2mat(Phi);

if size(Phi,2)>0.5*s*(s+1)
    s = ceil(0.5*(sqrt(1 + 8*size(Phi,2)) - 1)); Phi = {};
    for k = 1:p
        Phi{k} = cell2mat(rootTree(k,s));
    end
end
if iscell(Phi)
    Phi = cell2mat(Phi);
end
gamma = Phi(end,:)'; Phi = Phi(1:(end-1),:)';

if ~exist('method','var')
    method = 'L';
end

A = sym('a%d%d',[s s]); A = tril(A,-1); 
b = sym('b%d',[s 1]); 
c = sym('c%d',[s 1]); c(1) = 0;
assume([symvar(A) symvar(b) symvar(c)],'real');
assumeAlso([symvar(b) symvar(c)],'positive');
numIndeVars = 0.5*s*(s+1)-size(Phi,2);

if any(strcmpi(method,{'Gauss','G'}))
    
    
elseif any(strcmpi(method,{'Radau','R'}))
    

elseif any(strcmpi(method,{'Lobatto','L'}))
    

elseif any(strcmpi(method,{'Custom','C'}))
    if length(alpha)~=numIndeVars
        error(['Order conditions require ',numstr(numIndeVars),' free variables.']);
    end

end

depVars = [symvar(A) symvar(b) symvar(c)];
S = solve([Phi*b==gamma; c==A*ones(s,1)],depVars,'ReturnConditions',true,'Real',true);
depVars_str = string(depVars); indeVars = S.parameters; indeVars_str = string(indeVars);