function [ accept,accept1,accept2,stats ] = necSufCnd( S,A,gamma,Gamma )
% [ accept ] = percBasedSufCnd( S,A,gamma,Gamma )
%
% This function measures the sufficient condition
%   Agcmp'*Pg_ort*s < gamma*(1-sum(Agcmp'*pinv(Ag)',2))
%
% AND
%
% inv(Ag'*Ag) * (Ag'*s-gamma) > 0
%
% where 
% A: library matrix, 
% Ag is a submatrix of A with column indices correspond to Gamma, 
% Agcmp is a complement of Ag,
% Pg_ort: orthogonal projector onto the complement of Ag, 
% s: input signal, 
% gamma: a trade-off parameter, and
%
%
% Inputs
%   S : input signals [L,N] each column vectors are input signals
%   A : data matrix [L,p] each column vectors ar atoms
%   Gamma : set of indices (boolean with size [1,p] or integers 
%           in the range[1,p])
% Outputs
%   accept : boolean array [1,p]: True->the condition met False->not

[L,N] = size(A);
if islogical(Gamma)
    Gamma = find(Gamma);
    Gamma = reshape(Gamma,1,length(Gamma));
end
cGamma = setdiff(1:N,Gamma);
% 
Ag = A(:,Gamma);
Agcmp = A(:,cGamma);
AgtAg = Ag'*Ag;
% AgtAginv = inv(AgtAg);
Agpinv = pinv(Ag);
Pg_ort = eye(L) - Ag*Agpinv;

% Agcmp'*Pg_ort*s < gamma*(1-sum(Agcmp'*pinv(Ag)',2))
lhs1 = Agcmp' * Pg_ort * S;
rhs1 = gamma * ( 1 - sum(Agcmp' * Agpinv', 2) );
accept1 = all(bsxfun(@lt, lhs1, rhs1), 1);

% inv(Ag'*Ag) * (Ag'*s-gamma) > 0
lhs2 = AgtAg \ (Ag'*S-gamma);
accept2 = all(lhs2>0, 1);

accept = and(accept1,accept2);
stats = [];
stats.elementwise1 = false(N,size(S,2));
stats.elementwise1(cGamma,:) = bsxfun(@lt, lhs1, rhs1);
stats.PSCj = zeros(N,1);
stats.PSCj(cGamma) = ( 1 - sum(Agcmp' * Agpinv', 2) );
stats.PSCj(Gamma) = nan;
stats.nonlinearityj = zeros(N,size(S,2));
stats.nonlinearityj(cGamma,:) = lhs1;
stats.nonlinearity = Pg_ort * S;
stats.minCoeffs = lhs2;
end

