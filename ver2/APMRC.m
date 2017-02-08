function [ apmrc,mcc,nscc,opt_stats ] = APMRC( Y,A,gamma,Gamma )
% [ accept ] = percBasedSufCnd( S,A,gamma,Gamma )
%
% This function tests the approximately perfect model recoverty condition 
% (APMRC) of the paratially Lagrangian form of non-negative lasso 
% (PL-Nlasso):-
%
%      minimize         1/2||y-Ax||^2 + gamma*1^T*x
%      subject to       x >= 0
%
% where y is the Lx1 observation vector, A is the LxN dictionary matrix, 
% x is the Nx1 abundance vector, and gamma is the trade-off parameter. 
% The APMRC is composed of the two conditions, minimum coefficient
% condition(MCC) and non-negative versus subset coherence condition (NSCC).
%
% 1) MCC:   pinv(AG)*y > gamma * inv(AG'*AG) 
% 2) NSCC:  AGcmp'*PG_ort*y < gamma*(1-sum(AGcmp'*pinv(AG)',2))
%
% where AG is the subdictioanry matrix whose columens are associated with
% the indices in the subset "Gamma" of the columns of A, PG_ort is the 
% orthogonal projector onto the orthogonal complement of the span fo AG,
% and AGcmp is the matrix whose columns are indexed by the complement set
% of "Gamma".
%
%  Parameters
%     Y : input signals [L,M] each column vectors are input signals. L
%         is the number of dimensions and M is the number of samples.
%     A : the dictionary matrix [L,N] each column vectors ar atoms. N is
%         the number of atoms in the dictionary.
%     gamma : scalar, the trade-off parameter for PL-NLasso
%     Gamma : the set of indices (boolean with size [1,p] or integers 
%           in the range[1,p])
%  Returns
%     apmrc : boolean array [1,M]: True->the conditions are met False->not
%     mcc   : boolean array [1 M]: show whether MCC holds for each point.
%     nscc  : boolean array [1 M]: show whether NSCC holds for each point.
%
%  Optional returns
%     opt_stats : a struct containing more detail values inside the APMRC
%         FIELDs
%           Gamma:           integer indices of the subset "Gamma"
%           cGamma:          The complement of "Gamma"
%           gamma:           the trade-off parameter
%                  - for MCC
%           XG:              [p,M] the least square solution of Y=AG*XG
%           minCoeffs:       gamma*inv(AG'*AG)
%           mccfull:         [p,M] The elementwise outcome of MCC
%                  - for NSCC
%           PSCj:            [Nc,1] Positive Subset Coherence PSC(Gamma,j)
%                            Nc is the number of elements in cGamma.
%           errorVec:        [L,M], P_ort * Y
%           pnonliniearityj: [Nc,M], AGcmp'*PG_ort*Y.
%           nsccfull:        [Nc,M], the elementwise outcome of NSCC
%

%% check the validity of the input parameters
if size(Y,1) ~= size(A,1)
    error('The dimension of the input signal and the dictionary do not match\n');
end
[L,N] = size(A);

if ~isnumeric(gamma)
    error('The Value of the gamma is invalid.\n');
else
    if gamma<0
        error('The trade-off parameter gamma needs to be non-negative\n');
    end
end

if islogical(Gamma)
    Gamma = find(Gamma);
else
    if size(Gamma,1)>1 && size(Gamma,2)>1
        error('The shape of Gamma needs to be 1-dimension array\n.');
    elseif any(Gamma==0)
        error('Probably you forget to change "Gamma" to a logical array\n');
    end
end
Gamma = reshape(Gamma,1,length(Gamma));
p = size(Gamma,2);

if any(Gamma<1) || any(Gamma)>N
    error('The subset "Gamma" is invalid.\n.');
end

if length(Gamma) > L
    error('The size of "Gamma" is too large.The subdictionary is no more linear independent.\n');
end

%% compute several metrics.
cGamma = setdiff(1:N,Gamma);
AGcmp = A(:,cGamma);
AG = A(:,Gamma);
if rank(AG)<p
    error('The subdictionary is not full-rank.\n');
end
AGtAG = AG'*AG;


% 1) MCC: pinv(AG)*y > gamma * inv(AG'*AG)
% coefficient vector of the linear approximation of y on the span of AG
XG = AG\Y; 
mc = AGtAG\ones([p,1]) * gamma; % minimum coefficient
% cG = AGtAG \ (AG'*Y-gamma);
mccfull = bsxfun(@gt,XG,mc); % mcc conditions [p,M]
mcc = all(mccfull,1);

% 2) NSCC: AGcmp'*PG_ort*y < gamma*(1-sum(AGcmp'*pinv(AG)',2))
% orthogonal projection of Y onto the complement of the span of AG.
PG_ortY = Y - AG*XG;
% the left handside of NSCC ([cGamma,M])
pnonlinearity = AGcmp' * PG_ortY;
% Positive subset coherence
PSCs = (1 - sum(AG\AGcmp, 1))'; % [cGamma,1]
gPSCs = gamma * PSCs;
nsccfull = bsxfun(@lt, pnonlinearity, gPSCs); % [cGamma,M]
nscc = all(nsccfull,1);


apmrc = and(mcc,nscc);
opt_stats = [];
opt_stats.Gamma = Gamma;
opt_stats.cGamma = cGamma;
opt_stats.gamma = gamma;
% for MCC
opt_stats.XG = XG;
opt_stats.minCoeffs = mc;
opt_stats.mccfull = mccfull;
% for NSCC
opt_stats.PSCj = PSCs;
opt_stats.errorVec = PG_ortY;
opt_stats.pnonlinearityj = pnonlinearity;
opt_stats.nsccfull = nsccfull;
end

