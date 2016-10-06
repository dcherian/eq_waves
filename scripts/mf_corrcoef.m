function [r,rsig] = mf_corrcoef(x,y)
% FUNCTION MF_CORRCOEF Correlation and minimum correlation for 95% significance
% function [r,rsig] = mf_corrcoef(x,y)
% calculates results of corrcoef(x,y)
% using the non-NaN values of x,y
% and returns rsig, the minimum r that is significant at the 95% confidence level.
%
% Melanie Fewings, August 2008
% uses code from Carlos Moffat (rsignificant.m)
%
% INPUTS:
% x, y = vectors of data (all data points should be independent)
%
% OUTPUT:
% r = correlation coefficient from corrcoef.m
% rsig = minimum correlation for 95% significance

ii = find(~isnan(x+y));
% disp(['mf_corrcoef using ' num2str(length(ii)) ' of ' num2str(length(x)) ' points'])

if length(ii)>=3
    r = corrcoef(x(ii),y(ii));
    r = r(2,1);

    rsig = mf_rsignificant(length(ii));
else
    r = NaN;
    rsig = NaN;
end