function [ percamax_mrc,percamax_cnd,mcc ] = PERCMAX_MRC( Y,A,gamma,Gamma )
% [ percamax_mrc,percamax_cnd,mcc ] = PERCMAX_MRC( Y,A,gamma,Gamma )
%
% This function tests the PERC-MAX based model recovery condition (MRC) of
% the partially Lagrangian form of non-negative lasso (PL-Nlasso):-
%
%      minimize         1/2||y-Ax||^2 + gamma*1^T*x
%      subject to       x >= 0
%
% where y is the Lx1 observation vector, A is the LxN dictionary matrix, 
% x is the Nx1 abundance vector, and gamma is the trade-off parameter. 
% The PERC-MAX MRC is composed of the two conditions, minimum coefficient
% condition(MCC) and PERC-AMAX condition.
%
%  1) PERC-MAX condition: max(abs(A'*PG_ort*y),[],1) < gamma*PERC(Gamma)
%  2) MCC: inv(AG'*AG) * (AG'*y-gamma) > 0
%
% where AG is the subdictioanry matrix whose columens are associated with
% the indices in the subset "Gamma" of the columns of A, PG_ort is the 
% orthogonal projector onto the orthogonal complement of the span fo AG and
% PERC is short for positive exact recovery coefficient defined as:
% 
%     PERC(Gamma) =      min      1-sum(pinv(AG)*Aj)
%                    j \in Gcmp  |__________________|
%                                         ||
%                                         \/
%                                     PSC(Gamma;j)
%
% where AGcmp is the matrix whose columns are indexed by the complement set
% of "Gamma" and PSC is short for positive subset coherence.
%
% Parameters
%   Y : input signals [L,M] each column vectors are input signals
%   A : data matrix [L,N] each column vectors are atoms
%   gamma : trade-off parameter
%   Gamma : set of indices (boolean with size [1,N] or integers 
%           in the range[1,N])
% Outputs
%   percmax_mrc : boolean array [1,M], evaluating PERC-MAX & MCC
%                  True->the condition met False->not
%   percmax_cnd : boolean array [1,M], evaluating PERC-MAX condition
%                  True->the condition met False->not.
%   mcc          : boolean array [1,M], evaluating MCC
%                  True->the condition met False->not,
%
%   REFERENCE
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
% AgtAginv = inv(AgtAg);
Agpinv = pinv(Ag);
Pg_ort = eye(L) - Ag*Agpinv;

% ||A'*Pg_ort*s||_{infty} < gamma*PERC(Gamma)
lhs1 = max(abs(A'*Pg_ort*Y),[],1);
perc = PERC(A,Gamma);
percamax_cnd = (lhs1 < gamma*perc);

% inv(Ag'*Ag) * (Ag'*s-gamma) > 0
lhs2 = AgtAg \ (Ag'*Y-gamma);
mcc = all(lhs2>0, 1);

percamax_mrc = and(percamax_cnd,mcc);

end

