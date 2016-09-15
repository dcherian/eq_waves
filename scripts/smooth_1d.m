%smooth_1d(data,N,'power/amp','window') smooths a data vector to yield an N-point 
%		half-power or half-amplitude cutoff period.
%
%	Usage:
%	[smoothed_data]=smooth_1d(data,N,'power/amp','window');
%	where 
%		data = vector of data to be smoothed
%		N    = number of samples for desired half-power or half-amplitude period (non-integer OK)
%		'power/amp' (optional) can be either 'power' or 'amp' (for half-power or half-amplitude)
%			(default is 'amp')
%		'window' (optional) can be 'rect' (rectangle),'parzen','gauss', or 'cos2' (Hanning squared cosine)
%			(default is 'parzen')
%
%
%       2015, Tom Farrar, jfarrar@whoi.edu
%
% This function was modeled after run_avg.m, smoother_responses_v2.m, and Harris (1978, the window paper).
% Note that it can be used as a drop-in replacement for run_avg if the N-argument given to run_avg is divided by 24 for smooth_1d.


function [swav]=smooth_1d(srad,N,varargin)

if length(varargin)==0
  cut='amp';
  win='parzen';
elseif length(varargin)==1
  cut=varargin(1);
  win='parzen';
elseif length(varargin)==2
  cut=varargin(1);
  win=varargin(2);
end  

%desired cutoff freq:
fc=1./N;

%This is not the fastest way to set up the windows, but it is the most direct way I know of setting up the windows with the desired cutoff frequency.
% This is based off of smoother_responses_v2.m, and uses -3dB and -6dB points from Harris (1978) window paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First set up a vector of distance from estimation point: (in grid points)
dx=1;%0.01;
L=900;%30;
x=-L:dx:L;

if strcmp(win,'cos2')
%Cosine^2 filter
if strcmp(cut,'power')
Tcos2=1.44./fc/2;%-3 dB
elseif strcmp(cut,'amp')
Tcos2=2./fc/2;%-6 dB
end%amp/power
wcos2=cos(x.*pi./Tcos2).^2;
ff=find(-Tcos2/2>x|x>Tcos2/2);
wcos2(ff)=[];
coh_gain=0.5;
w=wcos2./(Tcos2*coh_gain);

elseif strcmp(win,'gauss')
%Gaussian window:
a=3;
if strcmp(cut,'power')
Tgauss=1.55./fc/2;%-3 dB
elseif strcmp(cut,'amp')
Tgauss=2.26./fc/2;%-6 dB; Harris says 2.18, which gives 52.5% power; 
end%amp/power
wgauss=exp(-0.5.*(a.*x./(Tgauss./2)).^2);
% I don't really understand how to truncate the window, I guess Tgauss/2
ff=find(-Tgauss/2>x|x>Tgauss/2);
wgauss(ff)=[];
coh_gain=0.4167;
w=wgauss./(Tgauss*coh_gain);

elseif strcmp(win,'parzen')
%Parzen window:
if strcmp(cut,'power')
Tparzen=1.82./fc/2;%-3 dB
elseif strcmp(cut,'amp')
Tparzen=2.55./fc/2;%-6 dB
end%amp/power
wparzen=0.*x;
ffi=find(0<=abs(x)&abs(x)<=Tparzen./4);
arg=abs(x)./(Tparzen/2);
wparzen(ffi)=1-6.*(arg(ffi).^2)+6.*(arg(ffi).^3);
ffi=find(Tparzen./4<=abs(x)&abs(x)<=Tparzen./2);
wparzen(ffi)=2.*(1-arg(ffi)).^3;
ff=find(wparzen==0);
wparzen(ff)=[];
coh_gain=0.3750;
w=wparzen./(Tparzen*coh_gain);

elseif strcmp(win,'rect')
%Rectangular wind0w:
if strcmp(cut,'power')
T1rm=0.89/fc/2;%-3dB
elseif strcmp(cut,'amp')
T1rm=1.21/fc/2;%-6 dB, Harris says 1.2, which gives 52% instead of 50% power; 1,21 gives 49.1%; 1.20001 gives 49.2%
end%amp/power
w1rm=0.*x+1;
ff=find(-T1rm/2>x|x>T1rm/2);%ff=find(-T1rm/2>=x|x>T1rm/2) exactly reproduces run_avg.m
w1rm(ff)=[];
coh_gain=1.000;
w=w1rm./(T1rm*coh_gain);
%(T1rm*coh_gain)

else error('Window specification not recognized')
end



%win=ones(M,1)./M;
[m,n]=size(srad);
if m==1 & n>1
  swav=conv2(1,w,srad,'same');
elseif n==1 & m>1
  swav=conv2(w,1,srad,'same');
else error('Input data is not a vector... aborting')
  return
end


