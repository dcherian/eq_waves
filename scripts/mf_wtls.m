function rr = mf_wtls(x,y,xerr,yerr,display_output)
% MF_WTLS Weighted total least squares line fit with NaNs
%
% function rr = mf_wtls_regress(x,y,xerr,yerr,display_output)
% calculates slope, intercept, and correlation coefficient
% rr = [slope,slope_err,int,int_err,r,rsig]
% using wtls_line.m
%
% Melanie Fewings, August 2008
% calls wtls_line.m
% (from http://www.mathworks.com/matlabcentral/fileexchange/17466)
% and uses code from nlinregi.m (Steve Lentz) to find 95% confidence
% intervals.
%
% INPUTS:
% x, y = two vectors of data to regress (find best line fit to y vs. x)
% xerr, yerr = uncertainties in x, y (can be a single value or a vector same size as x)
% NOTE: all data points in x, y should be independent but can contain NaNs
%
% OPTIONAL INPUTS:
% display_output = 1 for plot, 0 for no plot (default 0)
%
% OUTPUTS:
% rr = [slope,slope_err,int,int_err,r,rsig]
% slope, int = slope and intercept of the best-fit line
% slope_err, int_err = 95% confidence intervals for the slope and intercept
% r = linear correlation coefficient of x and y
% rsig = minimum correlation that is significant at 95% level with the
%        available number of data points

if nargin<5
    display_output = 0;
end

if sum(size(xerr))==2
    xerr = xerr.*ones(size(x));
elseif length(xerr)~=length(x)
    error('mf_wtls error: length(xerr) must equal 1 or length(x)')
end
if sum(size(yerr))==2
    yerr = yerr.*ones(size(x));
elseif length(yerr)~=length(y)
    error('mf_wtls error: length(yerr) must equal 1 or length(y)')
end

% [r,rsig] = mf_corrcoef(x,y);

ii = find(~isnan(x+y));
N = length(ii);

[slp,int,alpha,p1,chiopt,Cab,Calphap] = wtls_line(x(ii),y(ii), ...
                                                  xerr(ii),yerr(ii),0);
stdslp = sqrt(Cab(1));
stdint = sqrt(Cab(2));
% taken from nlinregi.m by Steve Lentz:
% now compute 95% ci's
st=tinv(0.975,N-2);
slperr=st*stdslp;
interr=st*stdint;

if display_output == 1
    figure; hold on;
    plot(x,y,'.','MarkerSize',20,'Color',.5.*[1 1 1])
    xlabel('X')
    ylabel('Y')
    ax = axis;
    title(['Weighted total least squares fit of Y on X with maximum 95% CIs'])
    plot(xlim,xlim.*slp+int,'k--','LineWidth',2)
    plot([0 max(xlim)],[0 max(xlim)].*(slp+slperr)+int+interr,'k--','LineWidth',0.75)
    plot([0 max(xlim)],[0 max(xlim)].*(slp-slperr)+int-interr,'k--','LineWidth',0.75)
    plot([min(xlim) 0],[min(xlim) 0].*(slp+slperr)+int-interr,'k--','LineWidth',0.75)
    plot([min(xlim) 0],[min(xlim) 0].*(slp-slperr)+int+interr,'k--','LineWidth',0.75)
    axis(ax);
end

rr = [slp,slperr,int,interr]; %,r,rsig];
