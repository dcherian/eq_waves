%% typical mode shapes

woa = load('../data/woa05.mat');

lon = 180; lat = 2;
ilon = find_approx(woa.X, lon, 1);
ilat = find_approx(woa.Y, lat, 1);

Zmode = avg1(woa.Z); % Z from woa05.mat
T = woa.temp(:,ilat,ilon);
S = woa.sal(:,ilat,ilon);
if all(isnan(T)) | all(isnan(S))
    error('no climatological data!');
end
dtdz = diff(T)./diff(woa.Z);

N2 = bfrq(S,T,woa.Z, lat);%10^(-6)*ones(32,1);
[Vmode, Hmode, c] = vertmode(N2,woa.Z,3,0);
% calculate temperature mode shape
Tmode = Vmode .* repmat(dtdz,1,size(Vmode,2));

figure;
clf;
plot(bsxfun(@rdivide, Tmode, -max(abs(Tmode))), -1*Zmode);
linex(0);
ylim([-2000 0]);
ylabel('Z (m)');
title(['Baroclinic mode shapes at (' num2str(lon) 'E, ' ...
       num2str(lat) 'N)']);
legend('1', '2', '3', 'Location', 'SouthEast');
beautify;
set(gcf, 'Position', [675 23 749 872]);
ax = gca;
ytick = ax.YTick;
yticklab = ax.YTickLabels;
liney([-1 -25 -50 -75 -100 -125 -150 -200 -250 -300 -500 -750]);
ax.YTick = ytick;

export_fig images/baroclinic-mode-shapes.png

%% test band pass filtering

fnameh = {'../data/dynht/dyn5s170w_dy.cdf';
          '../data/dynht/dyn2n95w_dy.cdf';};
%fnameh = ['../data/temp/t5n180w_dy.cdf'];
%dht = double(addnan(squeeze(ncread(fnameh,'T_20')),100))';

%dht = dht(:,5);
filt.halfdef = 'power';
filt.cutoff = sort(2./[0.135 0.155]);
filt.debugflag = 0;

SegmentLength = 365;

figure; maximize;

for ii=1:2
    subplot(1,2,ii);

    lon = ncread(fnameh{ii}, 'lon');
    lat = ncread(fnameh{ii}, 'lat');
    dht = double(addnan(squeeze(ncread(fnameh{ii} ,'DYN_13')),1000))';

    hdht = PlotSpectrum(dht, SegmentLength);

    filt.window = 'rect';
    hrect = PlotSpectrum(BandPass(dht, filt), SegmentLength);

    filt.window = 'gauss';
    hgauss = PlotSpectrum(BandPass(dht, filt), SegmentLength);

    filt.window = 'parzen';
    hparzen = PlotSpectrum(BandPass(dht, filt), SegmentLength);

    filt.window = 'butterworth';
    hbutt = PlotSpectrum(BandPass(dht, filt), SegmentLength);
    hbutt.Color = 'k';

    %linex(1./filt.cutoff);
    legend('raw dyn ht', 'rect filt', 'gauss filt', 'parzen filt', ...
           'butterworth', 'Location', 'SouthWest');
    title(['dyn ht at (' num2str(lon) 'E, ' ...
           num2str(lat) 'N) | [' num2str(sort(filt.cutoff), '%.2f ') ']']);
    beautify;
    grid on;
    linex(0.15, [], 'k');
    linex([0.19 0.11 0.085], [], [1 1 1]*0.6);
end

%export_fig -r300 images/filt-compare-dyn-ht-butter.png

%% time series - unfiltered and filtered


fnameh = {'../data/dynht/dyn5s170w_dy.cdf';
          '../data/dynht/dyn2n95w_dy.cdf';};
%fnameh = ['../data/temp/t5n180w_dy.cdf'];
%dht = double(addnan(squeeze(ncread(fnameh,'T_20')),100))';

%dht = dht(:,5);

filt.window = 'butterworth';
filt.halfdef = 'power';
filt.cutoff = sort(2./[0.135 0.155]);
filt.debugflag = 0;

SegmentLength = 365;

figure; maximize;

for ii=1:2
    subplot(2,1,ii); hold on;

    lon = ncread(fnameh{ii}, 'lon');
    lat = ncread(fnameh{ii}, 'lat');
    dht = double(addnan(squeeze(ncread(fnameh{ii} ,'DYN_13')),1000))';

    plot(dht-nanmean(dht));
    plot(BandPass(dht, filt));
    legend('raw dyn ht', filt.window, 'Location', 'SouthWest');
    title(['dyn ht at (' num2str(lon) 'E, ' ...
           num2str(lat) 'N) | [' num2str(sort(filt.cutoff), '%.2f ') ']']);
    beautify;
end

%% idealized test of band pass filtering
% dynamic height ω-k spectrum shows peaks at 0.15 cpd (bc2m1) and
% 0.2 cpd. Let's create synthetic time series and test filtering
% window

filt.halfdef = 'power';
filt.cutoff = 2./[0.135 0.155];
filt.debugflag = 0;
nsmooth = 5;
tvec = [1:8000];
ts = sin(2*pi*0.08*tvec) ...
     + sin(2*pi*0.1*tvec) ...
     + sin(2*pi*0.15 * tvec) ...
     + sin(2*pi*0.2 * tvec) ...
     + rand(size(tvec)) + 1;

figure;
hts = PlotSpectrum(ts);

filt.window = 'rect';
hrect = PlotSpectrum(BandPass(ts, filt));

filt.window = 'gauss';
hgauss = PlotSpectrum(BandPass(ts, filt));

filt.window = 'butterworth';
hbutter = PlotSpectrum(BandPass(ts, filt));

uistack(hbutter, 'bottom');
% filt.window = 'parzen';
% hparzen = PlotSpectrum(BandPass(ts, filt));
% hparzen.YData = smooth(hparzen.YData, nsmooth);

%linex(1./filt.cutoff);
legend('synthetic series', 'rect filt', ...
       'gauss filt', 'butterworth', 'parzen filt', ...
       'Location','NorthWest');
title([num2str(nsmooth) ' point smoothed | [' ...
       num2str(filt.cutoff) ']']);
beautify;
grid on;
%linex([0.15 0.2]);

%% white noise idealized filtering

tvec = [1:8000];
ts = rand(size(tvec));
winds = {'rect', 'gauss', 'parzen', 'butterworth'};

hi = 2./0.135;
lo = 2./0.155;

filt.halfdef = 'power';
filt.cutoff = sort(2./[0.135 0.155]);
filt.debugflag = 0;

figure;
hax = packfig(2,2);
for ii=1:length(winds)
    filt.window = winds{ii};

    axes(hax(ii));
    PlotSpectrum(ts);
    filt.N = lo; PlotSpectrum(FilterSeries(ts, filt));
    filt.N = hi; PlotSpectrum(FilterSeries(ts, filt));
    %PlotSpectrum(smooth_1d(ts, hi, 'power', filt.window) ...
    %             - smooth_1d(ts, lo, 'power', filt.window));
    PlotSpectrum(BandPass(ts, filt));
    title([filt.window ' | [' num2str(hi, '%.2f') ' ' num2str(lo, '%.2f') ']']);
    linex([0.08 0.1 0.15 0.2]);

    grid on
end
linkaxes(hax, 'xy');
ylim([1e-10 1])

%% compare all inferred mode structures

load bc2m1.mat

figure; hold on;
for mm=1:length(modes.lon)
    for nn=1:length(modes.lat)
        plot(modes.InferredMode{mm,nn}, -modes.depth{mm,nn}, 'k.');
    end
end

set(gcf, 'Position', [516 121 492 674]);
set(gca, 'XAxisLocation', 'top');
ylabel('Z (m)');
xlabel('Mode amplitude');
xlim([-0.2 1.2]);
linex(0);
beautify;

export_fig -r300 images/all-inferred-modes.png

%% check temperature structures at given moorings

load bc2m1-butterworth.mat

lonidx = [11 8 11];
latidx = [4 4 1];

figure; hold on; maximize;
hax = packfig(1,2);
set(hax, 'XAxisLocation', 'top');

for ii=1:length(lonidx)
    label = sprintf('(%dE, %dN)', modes.lon(lonidx(ii)), ...
                    modes.lat(latidx(ii)));
    axes(hax(1)); hold on;
    h(ii) = plot(data.Twoa{lonidx(ii), latidx(ii)}, -data.Zwoa, '.-', ...
                 'MarkerSize', 16, 'DisplayName', label);

    axes(hax(2)); hold on;
    plot(gradient(data.Twoa{lonidx(ii), latidx(ii)}, -data.Zwoa), -data.Zwoa, '.-', ...
         'MarkerSize', 16);
end
linkaxes(hax, 'y');
ylim([-800 0]);
axes(hax(1));
ylabel('Depth (m)'); xlabel('Temp (C)');
legend('Location', 'SouthEast');
beautify;
axes(hax(2));
xlabel('dT/dz');
hax(2).XTick(1) = [];
beautify;

export_fig images/T-dTdz.png

%% TAO/TRITON array seasonal bias?

load bc2m1.mat

time = [];
for mm=1:size(data.dhtfilt, 1)
    for nn=1:size(data.dhtfilt, 2)
        tdht = datetime(data.timedht{mm,nn}, 'ConvertFrom', 'juliandate');
        dht = data.dhtfilt{mm,nn};
        time = [time; tdht(~isnan(dht))];
    end
end

figure;
histogram(month(time), 12);
title(['histogram of valid filtered dynamic height observations vs. ' ...
       'month'])
xlabel('Month number');
xlim([1 12]);

export_fig -r150 images/dynht-month-hist.png

%% Compare inferred mode at (170W, 8S) for different filtering

linewidth = 2;
mm = 6; nn = 7;

PlotMode('bc2m1-rect', mm, nn);
gauss = load('bc2m1-gauss.mat');
errorbar(gauss.modes.InferredMode{mm,nn}, ...
         gauss.modes.depth{mm,nn} * -1, ...
         gauss.modes.InferredModeError{mm,nn}, ...
         'horizontal', ...
         'Marker', '.', 'MarkerSize', 12, ...
         'LineStyle', 'none', 'LineWidth', linewidth, ...
         'DisplayName', 'gauss');

butter = load('bc2m1-butterworth.mat');
errorbar(butter.modes.InferredMode{mm,nn}, ...
         butter.modes.depth{mm,nn} * -1, ...
         butter.modes.InferredModeError{mm,nn}, ...
         'horizontal', ...
         'Marker', '.', 'MarkerSize', 12, ...
         'LineStyle', 'none', 'LineWidth', linewidth, ...
         'DisplayName', 'butterworth');

hleg = legend('Location', 'SouthEast');
hleg.String{1} = 'rect';
beautify;

export_fig images/mode-shape-170w-8s-rect-gauss-butter.png

%%

lonrange = 5:11; latrange = 3:5;

[opt, plotopt] = DefaultOptions;
opt.filter_temp = 1;
modes = InferModeShape(opt, lonrange, latrange);
PlotModeMap(plotopt, lonrange, latrange, modes, opt);

plotopt.plotWTLS = 0;
PlotModeMap(plotopt, lonrange, latrange, modes, opt);

ylim([-500 0]);
export_fig -r150 images/10-14-bc2m1-eq.png

