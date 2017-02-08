function [ accept,accept1,accept2 ] = PERCMAX_MRC( S,A,gamma,Gamma )
% [ accept ] = percBasedSufCnd( S,A,gamma,Gamma )
%
% This function measures the sufficient condition
%    max (A'*Pg_ort*s) < gamma*PERC(Gamma)
%
% AND
%
% inv(Ag'*Ag) * (Ag'*s-gamma) > 0
%
% where 
% A: library matrix, 
% Ag is a submatrix of A with column indices correspond to Gamma, 
% Pg_ort: orthogonal projector onto the complement of Ag, 
% s: input signal,
% gamma: a trade-off parameter, and
% PERC is the positive exact recovery coefficient.
%
%
% Inputs
%   S : input signals [L,N] each column vectors are input signals
%   A : data matrix [L,p] each column vectors ar atoms
%   gamma : trade-off parameter
%   Gamma : set of indices (boolean with size [1,p] or integers 
%           in the range[1,p])
% Outputs
%   accept : boolean array [1,p]: True->the condition met False->not

[L,~] = size(A);

% preprocessing
Ag = A(:,Gamma);
AgtAg = Ag'*Ag;
% AgtAginv = inv(AgtAg);
Agpinv = pinv(Ag);
Pg_ort = eye(L) - Ag*Agpinv;

% max(A'*Pg_ort*s) < gamma*PERC(Gamma)
lhs1 = max(A'*Pg_ort*S,[],1);
perc = PERC(A,Gamma);
accept1 = (lhs1 < gamma*perc);

% inv(Ag'*Ag) * (Ag'*s-gamma) > 0
lhs2 = AgtAg \ (Ag'*S-gamma);
accept2 = all(lhs2>0, 1);

accept = and(accept1,accept2);

end

