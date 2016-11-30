%% Example script to run sparse scale free cholesky estimation.
% 
% This script demonstrates how to run SFCHL code. 
% As a resut, a precision matrix is estimated from data, along with sparse Cholesky factors. 
% Several evaluation metrics (Jaccard index, precision, recall) are printed to stdout.
%
%% Parameters
% There are several parameters to be tweaked:
%     dLthreshold         - convergence threshold 
%     numIterations       - maximum number of iterations allowed
%     normalizeCovariance - flag indicating correlation matrix to be used insted of covariance (recommended)
%     columnReordering    - flag indicating wheter to re-label nodes using Approximate Mininum Degree algorithm (recommended)
%     lambda              - regularization parameter
%     poolSize            - size of Matlab parallel pool (not applicable if called from Octave)
%% Data
% Data used in this example is generated from multivariate normal distribution with sparse precision matrix U. 
% Precision matrix was previously generated using the preferential attachment process.
%
%% Output
%
%     L     - Cholesky factor
%     Ueval - estimated precision matrix
%
% Author: Vladisav Jelisavcic, Ivan Stojkovic
% Work address
% email: 
% Website: 
% July 2016; Last revision: 20-July-2016

% Check for compatibility.
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

% Model name.
modelName = 'SFCHL';

% Name of the .mat file containing the dataset.
datasetName = 'scalefree2K';

% Dataset path.
dataPath = ['../../data/' datasetName];


%% ----- PARAMETERS ----

% Sparseness hyper-parameter.
lambda  = 0.0775;

% Threshold for convergence (calculated as sum(dL.^2)). 
dLthreshold = 0.1.^6; 

% Max number of iterations.
numIterations = 100;

% Set true to normalize covariance matrix (to correlation matrix).
Params.normalizeCovariance = true;

% Set true to enable Minimal degree reordering preprocessing step. 
Params.columnReordering = true;

% Size of Matlab parpool.
Params.poolSize = 0;

%----- DATA ----

% Add dataset and common files to path
parentpath = cd(cd('..'));
addpath(parentpath, [parentpath '/SFCHL']);

% Load dataset
load([dataPath '.mat']);

% Use correlation matrix instead of covariance
if(Params.normalizeCovariance)
    % Normalize Y to unit sphere.
    Y = Y./repmat(sqrt(sum(Y.*Y,1)),size(Y,1),1);  
end;

%% ----- INIT ----

% Init Matlab parallel pool.
if (~isOctave && Params.poolSize > 0)  % Note: parfor doesn't work on Octave. 
    pool = parpool('local',Params.poolSize);
end;

% Cholesky matrix (parameter to be optimized).
L = speye(N); % initially identity

disp('Calculating active set...');
tic;

% Precalculated (multilevel) active set.
activeSet = covLambda(Y, lambda);

% Starting active set.
activeSet0 = activeSet{1};

% Time needed to find starting active set.
elapsedTimePreproc = toc;

disp('Starting optimization...');
tic

reorderedCols = 1:N;

% Column reordering.
if(Params.columnReordering)
    reorderedCols = symamd(activeSet0); 

    U = U(reorderedCols, reorderedCols);

    Y = Y(:,reorderedCols);
    
    for i=1:length(activeSet)
      activeSet{i} = activeSet{i}(reorderedCols, reorderedCols);
    end;
end;

elapsedTimeReordering = toc;
tic

%% ------ OPTIMIZATION ----

% Main optimization procedure.
[L,iterations, dL] = SFCHL(Y, L, lambda, numIterations, dLthreshold, activeSet); 

% Time spent doing coordinate descent.
elapsedTime = toc;

% Total time.
elapsedTotal = elapsedTime + elapsedTimePreproc + elapsedTimeReordering;

% Calculate precision matrix.
Ueval = L*L';

%% ----- EVALUATION ----

% Calculate evaluation metrics and print results. 
postprocess(modelName, datasetName, U, Ueval, lambda, elapsedTotal);

% Cleanup (clean path, parpool, etc.).
if (~isOctave && Params.poolSize > 0)
    delete(pool);
end;

rmpath(parentpath);
rmpath([parentpath '/SFCHL']);
