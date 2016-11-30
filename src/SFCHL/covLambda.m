function [activeSet, count] = covLambda (Y, lambda, batchSize)
%COVLAMBDA function for active set calculation.
%
% 
% Syntax:  [activeSet, count] = covLambda (Y, lambda)
%
% Inputs:
%    y              - m \times n data vector, where m is number of samples, and n is number of variables
%    lambda         - 1 \times k decreasing vector of thresholding parameters, 
%                       where lambda(1) is equal to regularization parameter, 
%                       and k is level count
%    batchSize      - optimal parameter specifing batch size for covariance calculation
%                       used for big data (should be small enough to fit in RAM); default is 2000.
%    
% Outputs:
%    activeSet      - sparse Cholesky factor 
%    iterations     - number of iterations executed
%
% Example: 
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author: Vladisav Jelisavcic
% Work address
% email: 
% Website: 
% July 2016; Last revision: 20-July-2016

    if ~exist('batchSize','var'), batchSize = 2000; end    % Batch size.

    % Check for compatibility.
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

    % Total number of threholding levels. 
    Nlevels = length(lambda);
    
    % Number of nodes.
    N = size(Y,2);

    % Number of samples.
    Nsamples = size(Y,1);

    % Active set.
	  activeSet = cell(1,Nlevels);
  
    for i=1:Nlevels
        activeSet{i} = sparse(N,N);
    end;
  
    % Initiate frequency of values counting
    bounds = -1:0.01:1;
    count = zeros(1,length(bounds)-1);
    
    for i=1:batchSize:N
    
        % Last batch is smaller.
        batch = min(N-i,batchSize);
        
        % Covariance batch.
        data = Y'*Y(:,i:i+batch);

      	% count frequency of values
        %count = count + histcounts(data,bounds); 
        
         % Calculate thresholded covariance for current batch.
        activeBatch = sparse(abs(data) > lambda(Nlevels)).*data;
        
        activeSet{Nlevels}(:,i:i+batch) = activeBatch;
        
        for t=2:Nlevels
            activeBatch = (abs(activeBatch)>lambda(Nlevels-t+1)).*activeBatch;
            
            activeSet{Nlevels-t+1}(:,i:i+batch) = activeBatch; 
        end;
        
        disp([num2str(toc) ' batch:' num2str(ceil(i/batchSize)) '/' num2str(ceil(N/batchSize))]);
        
        if (isOctave)
          fflush(stdout);  % Needed to flush output on Octave.
        end;
    end;
    
    

end
