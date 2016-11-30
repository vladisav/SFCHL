function [] = postprocess (modelName, datasetName, U, Ueval, lambda, elapsedTotal)
%POSTPROCESS auxilary function for displaying results.
%
% 
% Syntax:  [] = postprocess (modelName, datasetName, U, Ueval, lambda, elapsedTotal)
%
% Inputs:
%    modelName      - name of the optimization function used
%    datasetName    - name of dataset used
%    U              - ground truth precision matrix
%    Ueval          - estimated precision matrix
%    lambda         - regularization parameter used (including threshold values)
%    elapsedTotal   - total time used for processing
%    
% Outputs:
%    none
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
    if ~exist('append','var'), append = false; end

    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    
    disp('-------------------------');
    disp(['model: ' modelName ', lambda: ' num2str(lambda)]);
    disp(['dataset: ' datasetName]); 
    disp('-------------------------');

    % Take elements above main diagonal (1+)
    trueU = triu(U,1);
    evalU = triu(Ueval,1);

    % Calculate Jaccard Index on upper triangular elements only
    JI = full(sum(sum(evalU & trueU))/sum(sum(evalU | trueU)))

    % Calculate nnz on both estimated and true graphs.
    sparsityU = nnz(U);
    sparsityLL = nnz(Ueval);
    
    disp(['nnz(U) =  ' num2str(sparsityU)]); 
    disp(['nnz(LL) =  ' num2str(sparsityLL)]);
    
    positiveA = trueU~=0;
    positiveB = evalU~=0;

    N = size(U,1);
    
    TP = nnz(positiveA & positiveB);
    FN = nnz((positiveA - positiveB)>0);
    FP = nnz((positiveB - positiveA)>0);
    TN = N*(N-1)/2 - (TP+FN+FP);
    
    precision = TP / (TP + FP);
    recall = TP / (TP + FN);
    accuracy = (TP+TN)/(TP+TN+FP+FN);
    
    disp(['precision =  ' num2str(precision)]);
    disp(['recall =  ' num2str(recall)]);
    disp(['accuracy =  ' num2str(accuracy)]);

    disp(['Total time:  ' num2str(elapsedTotal)]);
    disp('-------------------------');

end
