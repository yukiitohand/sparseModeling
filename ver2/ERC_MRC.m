function [ erc_mrc,erc_cnd,erc_mcc ] = ERC_MRC( S,E,thetaGamma,A,gamma,Gamma )
% [ erc_mrc,erc_cnd,erc_mcc ] = ERC_MRC( S,E,thetaGamma,A,gamma,Gamma )
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
% Inputs
%   S : input signals [L,N] each column vectors are input signals
%   E : errors [L,N] in the model (s = A\theta + e)
%   thetaGamma : abundances [|Gamma|,N] in the model (s = A\theta + e)
%   A : data matrix [L,p] each column vectors ar atoms
%   gamma : trade-off parameter
%   Gamma : set of indices (boolean with size [1,p] or integers 
%           in the range[1,p])
% Outputs
%   erc_mrc : boolean array [1,N]: True->the condition met False->not
%   erc_cnd : boolean array [1,N]: True->the condition met False->not
%   erc_mcc : boolean array [1,N]: True->the condition met False->not.
%
%   REFERENCE
%     Y Itoh, MF Duarte, M Parente, "Performance guarantees for sparse 
%     regression-based unmixing", 7th Workshop on Hyperspectral Image and 
%     Signal Processing: Evolution in Remote Sensing (WHISPERS), June 2015.
%
%     Y. Itoh, M. F. Duarte and M. Parente, "Perfect Recovery Conditions 
%     for Non-negative Sparse Modeling," in IEEE Transactions on Signal 
%     Processing, vol. 65, no. 1, pp. 69-80, 1 Jan.1, 2017.
%     doi: 10.1109/TSP.2016.2613067
%   
%   Copyright © 2019 Yuki Itoh

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
erc_cnd = (lhs1 <= gamma*erc);

% thetaGamma \succeq \gamma ||(Ag'Ag)^{-1}||_{infty,infty} - pinv(Ag)e
rhs2 = gamma*operatorNorm(AgtAginv,'inf','inf') - Agpinv*E;
erc_mcc = all(thetaGamma >= rhs2,1);

% take AND of the two
erc_mrc = and(erc_cnd,erc_mcc);

end

