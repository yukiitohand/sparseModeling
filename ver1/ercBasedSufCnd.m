function [ accept,accept1,accept2 ] = ercBasedSufCnd( S,E,thetaGamma,A,gamma,Gamma )
% [ accept ] = ercBasedSufCnd( S,E,thetaGamma,A,gamma,Gamma )
%
% This function measures the sufficient condition
%    ||A'*Pg_ort*s||_{infty} <= gamma*ERC(Gamma)
%
% AND
%
% thetaGamma >= \gamma*||(Ag'Ag)^{-1}||_{infty,infty} - pinv(Ag)*e
%
% where 
% A: library matrix, 
% Ag is a submatrix of A with column indices correspond to Gamma, 
% Pg_ort: orthogonal projector onto the complement of Ag, 
% s: input signal, e is an error vector,
% thetaGamma is an abundance vector,
% gamma: a trade-off parameter, and
% ERC is the exact recovery coefficient (ERC).
%
%
% Inputs
%   S : input signals [L,N] each column vectors are input signals
%   E : errors [L,N] in the model (s = A\theta + e)
%   thetaGamma : abundances [|Gamma|,N] in the model (s = A\theta + e)
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
AgtAginv = inv(AgtAg);
Agpinv = pinv(Ag);
Pg_ort = eye(L) - Ag*Agpinv;

%  ||A'*Pg_ort*s||_{infty} <= gamma*ERC(Gamma)
lhs1 = max(abs(A'*Pg_ort*S),[],1);
erc = ERC(A,Gamma);
accept1 = (lhs1 <= gamma*erc);

% thetaGamma \succeq \gamma ||(Ag'Ag)^{-1}||_{infty,infty} - pinv(Ag)e
rhs2 = gamma*operatorNorm(AgtAginv,'inf','inf') - Agpinv*E;
accept2 = all(thetaGamma >= rhs2,1);

% take AND of the two
accept = and(accept1,accept2);

end

