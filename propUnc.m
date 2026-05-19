function [dF, Fx] = propUnc(F,X,x,dx,C,c)
% propUnc(F,x,dx,X,C,c) evaluates the uncertainty and value of a function
% near a set of inputs with given uncertainties
% 
% Graham Holt, March 2026. Updated May 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% propUnc(F,x,dx,X)
% propUnc(___,C,c)
% dF = propUnc(___)
% [dF, Fx] = propUnc(___)
% 
%% Description
%   F - function relating the variables "X", expressed as a string
%   X - variables within "F", expressed as a row vector of strings
%   x - values of "X" to be considered, with each column corresponding 
%       to a different variable in "X"
%   dx - uncertainty in the values of "X", with each column 
%       corresponding to a different variable in X
%   C - symbols used for "c" in "F", expressed as a row vector of
%       strings
%   c - values for constants in "F", with element corresponding to a 
%       different constant in "C" 
%
% propUnc(F,x,dx,X) takes the function "F", with variables "X", and 
% returns the uncertainty of the function with variable values 
% "x" and uncertainties "dx"
%
% propUnc(___,C,c) clears up some notation by defining constants used in
% the function "F" which don't need different values
%
% [dF, Fx] = propUnc(___) includes the output of the function values for
% the given inputs

X = str2sym(X); F = str2sym(F);

if exist('c','var') && exist('C','var')
    C = str2sym(C);
    F = subs(F,C,c);
end

for i = 1:length(X)
    dFdX = diff(F,1,X(i));
    for j = 1:size(x,1)
        dF(j,i) = double(subs(dFdX,X,x(j,:)))*dx(j,i);
    end
end
for i = 1:size(x,1)
   Fx(i,:) = double(subs(F,X,x(i,:))); 
end
dF = double(sqrt(sum(dF.^2,2)));