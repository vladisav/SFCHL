function [L, iterations, dL] = SFCHLpar(y, L, lambda, numIterations, dLthreshold, activeSet, options)
%SFCHL main optimization function for fast scale free GMRF learning package 
% optimized for running on Matlab Parallel Computing Toolbox.
%
% 
% Syntax:  [L, iterations, dL] = SFCHL(y, L, lambda, numIterations, dLthreshold, activeSet, options)
%
% Inputs:
%    y              - m \times n data vector, where m is number of samples, and n is number of variables
%    L              - initial value
%    lambda         - 1 \times k decreasing vector of thresholding parameters, 
%                       where lambda(1) is equal to regularization parameter, 
%                       and k is level count 
%    numIterations - maximum number of iterations to execute
%    dLthreshold   - convergence threshold
%    activeSet     - 1 \times k cell array containing sparse covariance matrix
%                        thresholded by corresponding lambda
%    options       - optional struct containing additional options: 
%                        options.verbosity is used for different levels of displaying output to stdout,
%                        allowed values are {none,info,debug}, default is info
%                        options.display sets the display frequency, default is n/20
%
% Outputs:
%    L             - sparse Cholesky factor 
%    iterations    - number of iterations executed
%    dL            - gradient at final iteration
%
% Example: 
%    
%
% Other m-files required: assembleDataset.m
% Subfunctions: none
% MAT-files required: none
%
% See also: SFCHL

% Author: Vladisav Jelisavcic
% Work address
% email: 
% Website: 
% July 2016; Last revision: 20-July-2016


if ~exist('numIterations','var'), numIterations = 100; end
if ~exist('dLthreshold','var'), dLthreshold = 0.1.^6; end
if ~exist('L','var'), L = speye(N); end
if ~exist('lambda','var'), lambda = 0.9; end
if ~exist('activeSet','var'), activeSet = sparse(abs(data)>lambda); end
if ~exist('options','var'), options = struct; options.verbosity = 'info'; options.display = size(y,2)/20; end

% Check for compatibility.
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

% Number of nodes.
N = size(y,2);

% Number of data samples.
M = size(y,1);

% Derivatives.
dL = sparse(N,N);

% Number of iterations needed to converge.
iterations = zeros(1,N);

% Vector of decreasing thresholding values, 
% first one is equal to the regularization parameter. 
lambdaVec = lambda;

% Regularization parameter.
lambda = lambdaVec(1);

% Define verbosity level.
switch(options.verbosity)
    case 'none'  
        verbosity = 0;
    case 'info'  
        verbosity = 1;
    case 'debug' 
        verbosity = 2;
end;

if(verbosity > 0)
    disp('Starting main optimization function ...');
end;

initialActiveSet = activeSet{1};

% Start tracking time.
tic

% Iterate over each column.
parfor j=1:N
      % Column of the Choleksy factor to be optimized.
      X = L(:,j);

      % Index of the current thresholding value.
      thresholdIndex = 1;

      % Find active rows. 
      activeRows = find(initialActiveSet(:,j));
      
      % L is always triangular matrix.
      activeRows = activeRows(activeRows>=j);
      
      % Calculate covariances corresponding to active rows.
      S = assembleDataset(y,activeRows);
      
      % Main loop.     
      for (k=1:numIterations)

            % Counts number of additions/removals in current cd sweep.       
            changed = 0;
    
            % 
            Xold = X;
        
            % Iterate over each active element in column j.    
            for t=1:length(activeRows)
                i = activeRows(t);
                
                Sii = S(i,i);

                a = X'*S(:,i) - Sii*X(i);
                              
                if(i == j)
                    X(i) = (-a + sqrt(a^2 + 4*Sii))/(2*Sii); % Update diagonal elements.
                
                else
                    if(X(i) == 0)
                    
                        if(abs(a) < lambda) % Check the active set condition.
                            continue;      
                        end;
                        
                        % Increment the added/removed variables count.
                        changed = changed + 1;
                        
                        % Update off-diagonal elements.
                        X(i) = -(a - sign(a)*lambda)/Sii;                  
                    
                    else
                        X(i) = -(a + sign(X(i))*lambda)/Sii;
                        
                        % If sign is to be changed, stop at zero now, deal with it later.
                        if(sign(X(i)) ~= sign(Xold(i)))
                            changed = changed + 1;
                            X(i) = 0;
                        end;
                    end;
                end;
                
            end;
            
            dX = X-Xold;
            
            % Check convergence condition.
            if(sum(dX.^2) < dLthreshold && changed == 0)
               break;
            end;
             
            colSum = sum(abs(X));
        
            % Check if active set should be expanded.
            if(thresholdIndex < length(lambdaVec) && lambda < lambdaVec(thresholdIndex) * colSum)
                           
                Cthreshold = lambda / colSum;
                                
                for(thresholdIndex = thresholdIndex:length(lambdaVec))
                    if(Cthreshold >= lambdaVec(thresholdIndex))
                        break;
                    end;
                end;
                
                % Find active rows. 
                activeRows = find(activeSet{thresholdIndex}(:,j));
                
                % L is always triangular matrix.
                activeRows = activeRows(activeRows>=j);
                
                % Calculate covariances corresponding to active rows.
                S = assembleDataset(y,activeRows);
                
                if (verbosity > 1)
                    disp(['Col ' num2str(j) ' num active ' num2str(nnz(activeRows)) ' sum ' num2str(colSum) ' Cthreshold ' num2str(Cthreshold) ' threshold ' num2str(lambdaVec(thresholdIndex))]);         
                end;
            end;

      end;
      
      % Set output variables.
      L(:,j) = X;
      dL(:,j) = dX;
      iterations(j) = k;
      
      if(verbosity > 0 && mod(j,options.display) == 0)

          disp(['Finished column: ' num2str(j) '.']); 
      end;
      
      if (isOctave)
          fflush(stdout);  % Needed to flush output on Octave.
      end;
      
end;


end
