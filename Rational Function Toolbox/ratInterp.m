function [N,D,convFlag] = ratInterp(x,f,order,convInt)
% ratInterp(x,f) generates the simplest interpolating polynomial function
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% ratInterp(x,f)
% [N,D,~] = ratInterp(___,order)
% [N,D,~] = ratInterp(___,convInt)
% [N,D,~] = ratInterp(___)
% 
%% Description
% ratInterp(x,f) returns the simplest polynomial which evaluates to "f" at
% "x"
%
% ratInterp(___,order) specifies an order to use for a rational interpolant
%
% ratInterp(___,convInt) ensures that there are no poles between the left
% and right columns

if ~exist('convInt','var')
    convInt = [];
end

if exist('order','var')
    R = [vander2(x,order(1)+1), -f.*vander2(x,order(2)+1)];
    
    coeff = null(R)'; I = find(abs(coeff)>1e-12,1,'last');
    N = coeff(1:(order(1)+1))/coeff(I);
    N = N(find(abs(N)>1e-12,1):end);
    
    D = coeff((end-order(2)):end)/coeff(I);
    D = D(find(abs(D)>1e-12,1):end);
    
    p = polyRoots(D)';
    if imag(p)~=0
        return;
    end
    
    if isempty(convInt)
        return;
    end

    if all((convInt(:,1) < p) & (p > convInt(:,2)))
        return;
    else
        error('Convergence not achieved with provided order.');
    end
end

for k = 0:(length(x)-1)
    order = [length(x)-1-k, k];

    [N,D,convFlag] = ratInterp(x,f,convInt,order);
    if convFlag
        break;
    end
end