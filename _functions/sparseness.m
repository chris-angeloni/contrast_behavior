function S = sparseness(x);

% function S = sparseness(x);
%  computes sparseness of a vector. close to 1 means a pure
%  response to one bin, close to 0 is uniform response across bins.
% 
% from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3485085/pdf/zns10170.pdf

a = (sum(x)./length(x)).^2 ./ sum((x.^2)./length(x));
S = (1 - a) ./ (1 - 1/length(x));

S(isnan(S)) = 0;