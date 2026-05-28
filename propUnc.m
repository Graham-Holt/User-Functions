function [dF, Fx] = propUnc(F,X,x,dx,method,C,c)
% propUnc(F,x,dx,X) evaluates the uncertainty and value of a function
% near a set of inputs with given uncertainties
% 
% Graham Holt, March 2026. Updated May 2026
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
%        columns ('abs' method only), corresponding to a different variable 
%        in X
%   C  - symbols used for "c" in "F", expressed as a row vector of
%        strings
%   c  - values for constants in "F", with element corresponding to a 
%        different constant in "C" 
%
% propUnc(F,x,dx,X) takes the function "F", with variables "X", and 
% returns the uncertainty of the function with variable values 
% "x" and uncertainties "dx"
%
% propUnc(___,method) sets whether or not the output uses the standard,
% symmetric RMS method using partial derivatives (default), or computes the
% asymmetric upper and lower bounds of the values ('abs')
% 
% propUnc(___,C,c) clears up some notation by defining constants used in
% the function "F" which don't need different values
%
% [dF, Fx] = propUnc(___) includes the output of the function values for
% the given inputs

X = str2sym(X); F = str2sym(F);

% Pre-inserts constants where possible
if exist('c','var') && exist('C','var')
    C = str2sym(C);
    F = subs(F,C,c);
end

% Evaluates function as the given inputs
for k = 1:size(x,1)
    Fx(k,:) = double(subs(F,X,x(k,:)));
end

if strcmpi(method,'abs')
    % Checks evaluates at all possible uncertainies
    for k = 1:(2^size(x,2))
        if size(dx,2)==size(x,2)
            % Uses symmetric uncertainties
            s = 2*(dec2bin(k-1,size(x,2))-'0')-1;
            for j = 1:size(x,1)
                Fdx(j,k) = double(subs(F,X,x(j,:) + s.*dx(j,:)));
            end
        else
            % Uses asymmetric uncertainties
            s = (dec2bin(k-1,size(x,2))-'0') - 1 + (2:2:size(dx,2));
            for j = 1:size(x,1)
                Fdx(j,k) = double(subs(F,X,x(j,:)+dx(j,s)));
            end
        end
    end
    dF = double([min(Fdx,[],2) max(Fdx,[],2)] - Fx);
else
    % Uses partial derivatives to approximate uncertainties
    for k = 1:length(X)
        dFdX = diff(F,1,X(k));
        for j = 1:size(x,1)
            dF(j,k) = double(subs(dFdX,X,x(j,:)))*dx(j,k);
        end
    end
    dF = double(sqrt(sum(dF.^2,2)));
end