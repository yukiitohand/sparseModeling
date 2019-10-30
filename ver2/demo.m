%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create matrix of atoms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NA = 300; % number of each atom
L  = 200; % number of dimension of each atom.

% creating a library of atoms. 
A = randn(L,NA);
A = A ./ sqrt(sum(A.^2,1)); % normalization.

setU = 1:NA; % the set of all indicies associated with atoms.

setG = randsample(setU,5); % the set of indicies with which mixture will be created.
setGc = setdiff(setU,setG); % complement set of G.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate mixed samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ny = 100; % number of samples.
X = rand(length(setG),Ny); % generate coefficients randomly.
Y = A(:,setG)*X;

% noise
N = 0.003*randn(L,Ny);

% add noise
YN = Y+N;
% YN = Y;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically performing signal recovery using PL-Nlasso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = 0.01;
[X_hat] = sunsal(A,YN,'POSITIVITY','yes','VERBOSE','no','ADDONE','no', ...
    'lambda', gamma,'AL_ITERS',2000, 'TOL', 1e-8);
% You have to set strict threshold value tol=1e-8 or you need to use other
% more exact solvers using interior point method.

% Evaluate correct atoms are retrieved or not.
numer_recovry_result = and(all(X_hat(setG,:)>0),all(X_hat(setGc,:)<=0));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate recovery conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the model recovery conditions.
[ apmrc,mcc,nscc,opt_stats ] = APMRC(YN,A,gamma,setG);
[ percamax_mrc,percamax_cnd,mcc_percamax ] = PERCAMAX_MRC( YN,A,gamma,setG );
[ percmax_mrc,percmax_cnd,mcc_percmax ] = PERCMAX_MRC( YN,A,gamma,setG );
[ erc_mrc,erc_cnd,erc_mcc ] = ERC_MRC( YN,N,X,A,gamma,setG );

% Show results
fprintf('APMRC matches        %3.0f %% to the numerical recovery result.\n',100*mean(apmrc==numer_recovry_result));
fprintf('PERCAMAX MRC matches %3.0f %% to the numerical recovery result.\n',100*mean(percamax_mrc==numer_recovry_result));
fprintf('PERCMAX MRC matches  %3.0f %% to the numerical recovery result.\n',100*mean(percmax_mrc==numer_recovry_result));
fprintf('ERC MRC matches      %3.0f %% to the numerical recovery result.\n',100*mean(erc_mrc==numer_recovry_result));
