function zi = loess1d(xd,dd,xi,sx,check)
% loess filter estimate at xi of data [xd,dd] 
% zi = loess1d(xd,dd,xi,sx,check)
% xd is the coordinate
% dd is the data
% xi is the target coordinate
% sx is the loess filter width in units of xd
%
% If a 5th input (check) is given (it's value is not relevant) then
% the test for unsafe extrapolation in the case of a one-sided
% distribution of data within distance sx of the estimation 
% point is DISABLED
%
% John Wilkin (http://marine.rutgers.edu/dmcs/ms615/wilkin/#_Toc70145622)
%
%The quadratic loess smoother
%
%An example of weighted least squares fitting to multivariate data is the quadratic loess smoother.
%
%In one dimension, this can be used to smooth, filter or interpolate a time series of values that may or may not be at regular sampling intervals.
%
%An application in two or more dimensions could be to produce a gridded climatology from a set of data observed at irregularly spaced locations and times, such as a set of shipboard hydrographic observations of temperature, salinity, or other ocean properties.
%
% The flexibility to arbitrarily choose a set of weights means the least squares fitting approach can be tailored to meet a users’ notion of what constitutes a rational weighting scheme based on some a priori knowledge of the processes being observed and modeled.
%
%An example of this is the Climatology of the Australian Regional Seas (CARS):
%
%Ridgway, K., J. R. Dunn and J. L. Wilkin (2002), Ocean interpolation by 4-dimensional weighted least squares: Application to the waters around Australasia, Journal of Atmospheric and Oceanic Technology, 19, 1357-1375. (pdf)
%
%Length and time scales of variability
%
%The quadratic loess smoother requires an a priori choice of the length or time scale to apply in the selection of the smoothing weights. The smoother can be interpreted as a filter, since the linear weighting procedure is effectively implemented as a convolution of the weights with the data. It is more general in the sense that the data do not have to be at regular intervals, because the weights are computed simply as a function of the normalized distance function r.  It has been found, empirically, that the effective cutoff frequency of f­c the quadratic loess smoother when it is interpreted as a filter is fc ˜  L-1 where L is the half width (i.e. the normalization scale) used in the loess smoother. 
%
%If the loess smoother is to be used to deliberately remove certain scales of variability (i.e. as a filter), then selection of L is straightforward. However, if the objective is to use the smoother to do the best possible job of interpolating gaps in the data, then the smoothing scale should be adapted to the natural length or time scales of variability in the underlying physical process being observed.

 


% preallocate output
zi = NaN*ones(size(xi));

% hold the indicies of valid data points - this will be clipped so that we
% only determine a smoothed value for points that were valid unsmoothed

datalocations = 1:length(dd);

% check for NaNs in the data
nodata = find(isnan(dd));
if ~isempty(nodata)
  dd(nodata) = [];
  xd(nodata) = [];
  datalocations(nodata) = [];
end
if isempty(dd)
  zi=NaN.*xi';
  %bell;disp('There are no data');return
end

% loop over index of input coordinates
for j=1:length(xi)
  
  xj = xi(j);

  % centre coordinates on x
  x = xd - xj;
  
  % distance metric
  r = abs(x/sx);
  
  % only use data within r<=1
  Q = find(r<=1);
  r = r(Q);
  x = x(Q);
  d = dd(Q);

  % convert to vectors
  x = x(:); r = r(:); d = d(:);

  if (min([length(find(x>0)) length(find(x<0))]) > 3) | nargin>4
  
    % compute local weighted least squares quadratic estimate (loess)
    
    % form matrix A of data coordinates
    A = [ones(size(x)) x x.^2];
    
    % calculate weights w
    w = (1-r.^3).^3;
    
    % apply weights
    A = w(:,ones([1 3])).*A;
    rhs = w.*d;
    
    % c is the vector of coefficients of the quadratic fit:
    % fit = c(1) + c(2)*x + c(3)*x^2
    
    % solve least-squares fit of A*c = d
    c = A\rhs;
    
    % evaluate fitted value:
    % formally, zj = [1 xj xj^2]*c, but since we moved
    % origin to xj the answer is just c(1)
    zj = c(1);
    
    zi(j) = zj;

  end
  
end

coordinates
    A = [ones(size(x)) x x.^2];
    
    % calculate weights w
    w = (1-r.^3).^3;
    
    % apply weights
    A = w(:,ones([1 3])).*A;
    rhs = w.*d;
    
    % c is the vector of coefficients of the quadratic fit:
    % fit = c(1) + c(2)*x + c(3)*x^2
    
    % solve least-squares fit of A*c = d
    c = A\rhs;
    
    % evaluate fitte