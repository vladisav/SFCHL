function [data] = assembleDataset (y, activeCols, data)
%ASSEMBLEDATASET function for on-the-fly calculation of required portion of covariance matrix. 
%
% 
% Syntax:  [data] = assembleDataset (y, activeCols, data)
%
% Inputs:
%    y              - m \times n data vector, where m is number of samples, and n is number of variables
%    activeCols     - vector containing indices of active elements in a single row
%    data           - sparse matrix containing required portion of covariance matrix
% 
% Outputs:
%    data           - n \times n sparse matrix containing a required portion of covariance matrix
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
      % Number of nodes.
      N = size(y,2);

      % Number of data samples.
      M = size(y,1);

      if ~exist('data','var'), data = sparse(N,N); end

      % FIXME add data re-use

      t = y(:,activeCols);

      S = t'*t;
      
      data(activeCols,activeCols) = S;
      
end
