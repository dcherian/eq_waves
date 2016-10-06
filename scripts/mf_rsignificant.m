function rsig = mf_rsignificant(num);
% FUNCTION MF_RSIGNIFICANT Lowest significant correlation coefficient
% function rsig = mf_rsignificant(num)
% Calculates the lowest correlation, rsig, that is significant at the
% 95% confidence level with "num" number of independent points
%
% Melanie Fewings, August 2008
% uses code from Carlos Moffat
%
% INPUT:
% num = number of independent data points
%
% OUTPUT:
% rsig: correlation must be higher than rsig 
% to be significant at the 95% confidence level
%
% see Emery & Thomson, Data Analysis Methods in Physical Oceanography, p. 253

if num>3
    r = .001:.001:.999;
    y = sqrt(num-3)./2.*log((1+r)./(1-r));
    rsig = interp1(y,r,1.96);
else
    rsig = NaN;
end



