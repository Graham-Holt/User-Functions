function [dF, Fx] = propUnc(F,X,x,dx,C,c,method)
% propUnc(F,x,dx,X) evaluates the uncertainty and value of a function
% near a set of inputs with given uncertainties
% 
% Graham Holt, March 2025. Updated July 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% propUnc(F,x,dx,X)
% propUnc(___,method)
% propUnc(___,C,c)
% dF = propUnc(___)
% [dF, Fx] = propUnc(___)
% 
%% Description
%   F  - function relating the variables "X", expressed as a string
%   X  - variables within "F", expressed as a row vector of strings
%   x  - values of "X" to be considered, with each column corresponding to 
%        a different variable in "X"
%   dx - uncertainty in the values of "X", with each column, or pair of 
%        columns (w/ 'abs' method), corresponding to a different variable 
%        in X
%   C  - symbols used for "c" in "F", expressed as a row vector of
%        strings
%   c  - values for constants in "F", with element corresponding to a 
%        different constant in "C" 
%
% propUnc(F,x,dx,X) takes the function "F", with variables "X", and 
% returns the uncertainty of the function with variable values 
% "x" and uncertainties "dx".
%
% propUnc(___,C,c) clears up some notation by defining constants used in
% the function "F" which don't need different values. If not included, use
% empty string or cell array.
%
% propUnc(___,method) sets whether or not the output uses the standard,
% symmetric RMS method using partial derivatives (default), or computes the
% asymmetric upper and lower bounds of the values ('abs').
%
% [dF, Fx] = propUnc(___) includes the output of the function values for
% the given inputs.

% Converts strings into usable symbolic functions
X = str2sym(X); F = str2sym(F); symmetric = size(x,2)==size(dx,2);

% Pre-inserts constants where possible
if exist('c','var') && exist('C','var')
    C = str2sym(C);
    F = subs(F,C,c);
end

% Selects default method
if ~exist('method','var')
    if ~symmetric
        method = 'abs';
    else
        method = 'rms';
    end
end

% Evaluates function at the given inputs
for k = 1:size(x,1)
    Fx(k,:) = double(subs(F,X,x(k,:)));
end

if strcmpi(method,'abs')
% Evaluates at all possible uncertainies for comparision

% Puts symmetric uncertainies into asymmetric form
if symmetric
    collate = reshape([0; size(x,2)]+(1:size(x,2)),1,[]);
    dx = [-dx dx]; dx = dx(:,collate);
end

for k = 1:(2^size(x,2))    
    % Selects each combination of uncertainties to comparison
    selectUnc = (dec2bin(k-1,size(x,2))-'0') - 1 + (2:2:size(dx,2));
    for j = 1:size(x,1)
        Fdx(j,k) = double(subs(F,X,x(j,:)+dx(j,selectUnc)));
    end
end
dF = double([min(Fdx,[],2) max(Fdx,[],2)] - Fx);

else
% Uses partial derivatives to approximate effect of uncertainities

% Considers largest deviation for uncertainty
if ~symmetric
    dx = cat(3,dx(:,1:2:end),dx(:,2:2:end));
    dx = max(abs(dx),[],3);
end

for k = 1:length(X)
    for j = 1:size(x,1)
        dF(j,k) = subs(diff(F,1,X(k)),X,x(j,:))*dx(j,k);
    end
end
dF = double(sqrt(sum(dF.^2,2)));

end