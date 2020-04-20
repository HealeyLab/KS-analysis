function [output] = classifyKNN(obj, searchD, searchC, testD)
%CLASSIFYKNN Performs k-NN classification
%   Using training raw data (searchD), the categories that correspond to
%   the training data (searchC), and the test data (testD), classifies the
%   testD.

% search data (searchD), search data classifications (searchC) 
% and test data and spits out an array in which 
[NN, ~]=knnsearch(searchD, testD, 'K', 1);
% searchC(NN) spits out an m x 3 matrix of the classifications
output = searchC(NN);

end

