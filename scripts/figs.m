%% test band pass filtering

fnameh = ['../data/dynht/dyn5n180w_dy.cdf'];
%fnameh = ['../data/temp/t5n180w_dy.cdf'];
dht = double(addnan(squeeze(ncread(fnameh,'DYN_13')),1000))';
%dht = double(addnan(squeeze(ncread(fnameh,'T_20')),100))';

%dht = dht(:,5);
filt.halfdef = 'power';
filt.cutoff = sort(2./[0.14 0.15]);
filt.debugflag = 0;
nsmooth = 5;

figure;
hdht = PlotSpectrum(dht);
hdht.YData = smooth(hdht.YData, nsmooth);

filt.window = 'rect';
hrect = PlotSpectrum(BandPass(dht, filt));
hrect.YData = smooth(hrect.YData, nsmooth);

filt.window = 'gauss';
hgauss = PlotSpectrum(BandPass(dht, filt));
hgauss.YData = smooth(hgauss.YData, nsmooth);

filt.window = 'parzen';
hparzen = PlotSpectrum(BandPass(dht, filt));
hparzen.YData = smooth(hparzen.YData, nsmooth);

%linex(1./filt.cutoff);
legend('DHT', 'rect filt', 'gauss filt', 'parzen filt');
title([num2str(nsmooth) ' point smoothed dyn ht | [' num2str(1./filt.cutoff) ']']);
beautify;
grid on;
linex(0.15, 'bc2m1', 'k');
linex([0.19 0.11 0.085], [], [1 1 1]*0.5);

%% idealized test of band pass filtering
% dynamic height ω-k spectrum shows peaks at 0.15 cpd (bc2m1) and
% 0.2 cpd. Let's create synthetic time series and test filtering
% window

filt.halfdef = 'power';
filt.cutoff = 2./[0.14 0.15];
filt.debugflag = 0;
nsmooth = 5;
tvec = [1:8000];
ts = sin(2*pi*0.08*tvec) ...
     + sin(2*pi*0.1*tvec) ...
     + sin(2*pi*0.15 * tvec) ...
     + sin(2*pi*0.2 * tvec) ...
     + rand(size(tvec));

figure;
hts = PlotSpectrum(ts);

filt.window = 'rect';
hrect = PlotSpectrum(BandPass(ts, filt));
hrect.YData = smooth(hrect.YData, nsmooth);

filt.window = 'gauss';
hgauss = PlotSpectrum(BandPass(ts, filt));
hgauss.YData = smooth(hgauss.YData, nsmooth);

% filt.window = 'parzen';
% hparzen = PlotSpectrum(BandPass(ts, filt));
% hparzen.YData = smooth(hparzen.YData, nsmooth);

%linex(1./filt.cutoff);
legend('synthetic series', 'rect filt', ...
       'gauss filt', 'parzen filt', ...
       'Location','NorthWest');
title([num2str(nsmooth) ' point smoothed | [' ...
       num2str(filt.cutoff) ']']);
beautify;
grid on;
%linex([0.15 0.2]);

%% white noise idealized filtering

tvec = [1:8000];
ts = rand(size(tvec));
winds = {'rect', 'gauss', 'parzen', 'cos2'};

hi = 2./0.14;
lo = 2./0.15;

figure;
hax = packfig(2,2);
for ii=1:length(winds)
    window = winds{ii};

    axes(hax(ii));
    PlotSpectrum(ts);
    PlotSpectrum(smooth_1d(ts, hi, 'power', window));
    PlotSpectrum(smooth_1d(ts, lo, 'power', window));
    PlotSpectrum(smooth_1d(ts, hi, 'power', window) ...
                 - smooth_1d(ts, lo, 'power', window));
    title([window ' | [' num2str(hi, '%.2f') ' ' num2str(lo, '%.2f') ']']);
    linex([0.08 0.1 0.15 0.2]);

    grid on
end
linkaxes(hax, 'xy');
