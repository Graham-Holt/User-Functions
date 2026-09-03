function varargout = designControl(varargin)
% designControl(A,B,C,D,eig,method) computes gain matrices for linearized
% dynamical systems
% 
% Graham Holt, August 2026. Updated September 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% designControl(A,B,C,D,eig,method)
% K = designControl(___)
% [K,L] = designControl(___)
% [K,L,Kr] = designControl(___)
% 
%% Description
% designControl(A,B,C,D,eig,method) takes in the linearized system 
% (A,B,C,D) and the desired eigenvalues to place (eig). "method" is used to 
% determine if the control uses state feedback ('x') or output feedback 
% ('y'). This form returns the control gain matrix.
%
% [K,L] = designControl(___) will include the state estimation gain matrix. 
% This will return an empty matrix if output feedback is used.
%
% [K,L,Kr] = designControl(___) will include the control gain matrix 
% attached to a desired reference away from the origin of linearization. 
% The reference can be taken to be a target state ('-x') or a target output 
% ('-y'). If not specified in "method", the default is a reference state.

if nargin < 6
    error('Too few input arguments.');
elseif nargin > 6
    error('Too many input arguments.');
end

A = varargin{1};
B = varargin{2};
C = varargin{3};
D = varargin{4};

[n,m] = size(B);
[p,~] = size(C);

eig = varargin{5};
method = lower(char(varargin{6}));

if isscalar(method)
    method = [method, 'x'];
else
    method = method(1:2);
end

switch method(1)
case 'x'
    K = placeX(A,B,eig(1:n));
    L = placeX(A',C',eig((n+1):end)).';
case 'y'
    K = placeY(A,B,C,D,eig);
    L = [];

    Kinf = K(1:m,1:p) - K(1:m,(p+1):end)/K((m+1):end,(p+1):end)*K((m+1):end,1:p);
otherwise
    error('Improper method chosen.');
end

if nargout > 2
    
switch method
case 'xx'
    Kr = -(A - B*K)\B;
case 'xy'
    Kr = D - (C + D*K)/(A - B*K)*B;
case 'yx'
    Kr = -(A - B*Kinf/(eye(p) + D*Kinf)*C)\B*(eye(m) + Kinf/(eye(p) + D*Kinf)*D);
case 'yy'
    Kr = -(A - B*Kinf/(eye(p) + D*Kinf)*C)\B*(eye(m) + Kinf/(eye(p) + D*Kinf)*D);
    Kr = (eye(p) + D*Kinf)\(D + C*Kr);
otherwise
    error('Improper method chosen.');
end
Kr = pinv(Kr);

end

switch nargout
case {0,1}
    varargout = {K};
case 2
    varargout = {K,L};
case 3
    varargout = {K,L,Kr};
otherwise
    error('Too many output arguments.');
end

end

function K = placeX(A,B,eig)
% placeX(A,B,eig) places poles for state feedback design.
% 
% Graham Holt, August 2026. Updated September 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% placeX(A,B,eig)
% K = placeX(___)
% 
%% Description
% placeX(A,B,eig) numerically solves for a gain matrix to place closed-loop 
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
        K(k,j,:) = reshape(S.("k"+num2str(k)+num2str(j)),1,1,[]);
    end
end
K = double(K);
[~,I] = mink(reshape(pagenorm(K),1,[]),1);
K = K(:,:,I);

end

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

end