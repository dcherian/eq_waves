%%
mm = 7; nn = 4;
opt = DefaultOptions;
opt.TagainstDHT = 1;
opt.numMC = 2000;

InferOneLocation(mm,nn,opt);

%% Test monte carlo with peaks


%% linear regression monte carlo

xx = randn([5000 1]);
yy = 0*xx + 4*randn([5000 1]);
coeff = dcregress(xx, yy);
m0 = coeff(2)
noise_estimate = yy - m0*xx;
[NoiseAmp, NoiseSlope] = EstimateNoiseSpectrum(noise_estimate, [], 0)

for ii =1:5000
    noisevec = synthetic_timeseries_known_spectrum(length(xx), 1, ...
                                                   NoiseAmp, NoiseSlope);
    yy1 = m0*xx + noisevec;
    coeff = dcregress(xx, yy1, 0, 0, 0);
    c(ii) = coeff(1);
    m(ii) = coeff(2);
end
histogram(m); linex(m0);

%% test interpgaplength
opt.InterpGapLength = 0;
modes = InferModeShape(opt, tao, mm, nn);
PlotMode(modes, mm, nn, plotopt);

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

%% compare correlations with amplitudes

% load('bc2m1-butterworth');

figure; hold on;
for ii=1:length(modes.lon)
    for jj=1:length(modes.lat)
        plot(modes.InferredModeOLS{ii,jj}, modes.corr{ii,jj}, ...
        'k.');
        if any(modes.InferredModeOLS{ii,jj} > 1)
            disp([ii jj])
        end
    end
end
xlabel('Inferred mode amplitude');
ylabel('Corr. coeff');

%% check regression slopes
mm = 8 ; nn = 4;

tao = ReadTaoTriton(mm, nn);
[opt, plotopt] = DefaultOptions;
opt.numMC = 1e3;
modes = InferModeShape(opt, tao, mm, nn);
PlotMode(modes, mm, nn);

%% slope burger number
lat = 2.5; % degrees
f = 2* (2*pi/86400) * sind(lat);
alpha = 0.01; % bottom slope
N = 1e-2;

(alpha * N/f)^2

%% correlation (white noise) significance test

tic;
numMC = 5e3;
numel = 8000;

rred = nan([numMC 1]);
rwhite = nan([numMC 1]);
for ii=1:numMC
    x = rednoise([numel 1]);
    y = rednoise([numel 1]);

    rmat = corrcoef(x,y);
    rred(ii) = rmat(1,2);

    x = whitenoise([numel 1]);
    y = whitenoise([numel 1]);

    rmat = corrcoef(x,y);
    rwhite(ii) = rmat(1,2);
end

rsort = sort(abs(rwhite));
sigwhite = rsort(floor(0.975*numMC));
rsort = sort(abs(rred));
sigred = rsort(floor(0.975*numMC));
disp(['white: ' num2str(sigwhite, '%.4f')]);
disp(['red: ' num2str(sigred, '%.4f')]);
corr_sig(numel-2, 0.95)

clf
hr = histogram(rred, 20, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
hw = histogram(rwhite, 20, ...
               'FaceColor', 'w');

hthsig = linex([-1 1]*corr_sig(numel-2, 0.95));
hwsig = linex([-1 1]*sigwhite, [], 'k');
hrsig = linex([-1 1]*sigred, [], 'r');

legend([hw hr hwsig{1} hrsig{1} hthsig{1}], ...
       'White', 'Red', 'White 95', 'Red 95', 'Theory 95');

set(gca, 'XTickMode', 'auto');
ylabel('Counts');
xlabel('Correlation coefficient');
title(['Num MC iterations = ' num2str(numMC) ' | Time series length = ' num2str(numel)]);
beautify;
export_fig images/monte-carlo-corr-coeff-white-red.png
toc;

%%
opt = DefaultOptions;
in = synthetic_timeseries_known_spectrum(1e5, 1, 1e-6, -5)/1e5;
figure;
hax(1) = subplot(211);
plot(in);
hax(2) = subplot(212);
plot(BandPass(in, opt.filt));
linkaxes(hax, 'x');

%% Butterworth filter step response
opt = DefaultOptions;
figure; hold on;
for ii=1:4
    [b,a] = butter(ii, sort(2./opt.filt.cutoff/(1/2)), 'bandpass');
    [h,t] = stepz(b,a);
    hplt(ii) = plot(t,h+(ii-1)/6, 'DisplayName', num2str(ii));
end
legend('Location', 'SouthEast')
title('(1,2,3,4) order butterworth filter step response uzing stepz')
resizeImageForPub('portrait');
beautify([11 12 13])

export_fig images/butterworth-step-response.png

%% compare windows
plotopt.window = 'rect';
PlotModeMap(plotopt);
export_fig -r300 images/09-28-bc2m1-rect.png

plotopt.window = 'gauss';
PlotModeMap(plotopt);

export_fig -r300 images/09-28-bc2m1-gauss.png

ylim([-500 0]);
export_fig -r300 images/09-28-bc2m1-top500.png

%% bc2m1 with different filtering
clear opt

% main options
opt.debug = 0; % debugging spectrum plots
opt.filter_temp = 1; % filter temperature also?

opt.filt.halfdef = 'power'; % how is filt.N defined?
opt.filt.N = NaN; % will be set based on cutoff later.
opt.filt.debugflag = 0; % debugging spectrum plots in BandPass()
opt.name = 'bc2m1';
opt.filt.cutoff = 2./[0.135 0.155]; % (days) band pass filter windows

opt.filt.window = 'rect';
InferModeShape(opt);

opt.filt.window = 'gauss';
InferModeShape(opt);

opt.filt.window = 'butterworth';
InferModeShape(opt);

%% bc2m2

opt.name = 'bc2m2';
opt.filt.cutoff = 2./[0.19 0.23]; % (days) band pass filter windows
InferModeShape(opt);

plotopt.name = opt.name;
plotopt.nmode = [1 2]; % which theoretical mode am I looking for?
plotopt.plotcorr = 0;
plotopt.plotstd = 0;
PlotModeMap(plotopt);

export_fig -r300 images/09-24-bc2m2.png

%% monte carlo regression slope

len = 5000;

% unfiltered white noise
disp('White noise')
[mWhite, slWhite, rWhite, wWhite] = TestMC(0,0,0,len);
fitWhiter = fitdist(wWhite, 'normal')
rdof = 1./fitWhiter.sigma^2+3
fitWhitet = fitdist(mWhite, 'tlocationscale')
fitWhiteg = fitdist(mWhite, 'normal')
% this should give standard t-values for large degrees of freedom
% calc95(mWhite)
% tinv(0.025, len)

% band passed white noise
disp('Filtered white noise')
[mWhiteFilter, slWhiteFilter, rWhiteFilter, wWhiteFilter] = TestMC(0,1,0);
fitWhiteFilterr = fitdist(wWhiteFilter, 'normal')
rdof = 1./fitWhiteFilterr.sigma^2+3
fitWhiteFiltert = fitdist(mWhiteFilter, 'tlocationscale')
fitWhiteFilterr = fitdist(wWhiteFilter, 'normal')
fitWhiteFilterg = fitdist(mWhiteFilter', 'normal')

% red noise
disp('Red noise')
[mRed, slRed, rRed, wRed] = TestMC(-3,0);
fitRedr = fitdist(wRed, 'normal')
rdof = 1./fitRedr.sigma^2+3
fitRedt = fitdist(mRed, 'tlocationscale')
fitRedg = fitdist(mRed, 'normal')

% band passed red noise
disp('Filtered red noise:')
[mRedFilter, slRedFilter, rRedFilter, wRedFilter] = TestMC(-2, 1);
fitRedFilterr = fitdist(wRedFilter, 'normal')
rdof = 1./fitRedFilterr.sigma^2+3
fitRedFiltert = fitdist(mRedFilter, 'tlocationscale')
fitRedFilterg = fitdist(mRedFilter, 'normal')

%%
nbins = 40;

figure;
histogram(slWhiteFilter, nbins, 'FaceColor', [1 1 1]*0.25);
hold on
histogram(slWhite, nbins, 'FaceColor', [1 1 1]);

figure;
cla;
histogram(mWhiteFilter, nbins, 'FaceColor', [1 1 1]*0.25, ...
          'EdgeColor', 'none');
hold on;
histogram(mWhite, nbins, 'FaceColor', [1 1 1], ...
          'EdgeColor', 'none');
histogram(mRed, nbins, 'EdgeColor', 'none')
histogram(mRedFilter, nbins, 'EdgeColor', 'none')
legend('Filtered white noise', 'White noise')

%% white noise Filtering
opt = DefaultOptions;

xx = synthetic_timeseries_known_spectrum(len, 1, 1, -3);
yy = synthetic_timeseries_known_spectrum(len, 1, 1, 0);

figure;
subplot(121);
plot(xx, yy, '*');
[coeff,~,~,err] = dcregress(xx, yy, [], 0, 0, 0)

xx = BandPass(xx, opt.filt);
yy = BandPass(yy, opt.filt);
subplot(122);
plot(xx, yy, '*');
[coeff,~,~,err] = dcregress(xx, yy, [], 0, 0, 0)

%% correlation structure and filtering

figure('Position', [146 1030 560 573]);
hax(1) = subplot(311);
CorrFilter(0, hax(1));
xlabel('');
beautify([13 14 15]);

hax(2) = subplot(312);
CorrFilter(-3, hax(2));
xlabel('');
beautify([13 14 15]);

hax(3) = subplot(313);
CorrFilter(-5, hax(3));
beautify([13 14 15]);

linkaxes(hax, 'x');
xlim([-1 1]*30)

export_fig -transparent images/correlation-structure-filtering.png

%% figure out data frequency slopes

SubsetLength = 256;
if ~exist('tao', 'var'), tao = ReadTaoTriton; end

Tslope = nan([11 7]);
for mm=1:11
    for nn=1:7
        if isempty(tao.dht{mm,nn}), continue; end

        [Sdht, freq] = GappySpectrum(tao.dht{mm,nn}, SubsetLength);

        St = nan([length(freq) size(tao.T{mm,nn}, 1)]);
        for zz=1:size(tao.T{mm,nn}, 1)
            try
                St(:,zz) = GappySpectrum(tao.T{mm,nn}(zz,:), SubsetLength);
            catch ME
                St(:,zz) = NaN;
            end
        end

        % avg over depth
        St = nanmean(St, 2);

        [coeff, conf] = dcregress(log(freq), log(St), length(freq), [], 0);
        Tslope(mm,nn) = coeff(2);
        [coeff, conf] = dcregress(log(freq), log(Sdht), length(freq), [], 0);
        DHTslope(mm,nn) = coeff(2);
    end
end

figure;
histogram(DHTslope);
hold on;
histogram(Tslope);
legend('Dyn ht', 'temp');

%% monte carlo distribution fits

load montecarlo

figure;
hax(1) = subplot(121);
FitAndPlotDist(slfit, [], hax(1));
title('Regression slope');
pbaspect([1.3 1 1]);
beautify([12 14 15]);

hax(2) = subplot(122);
FitAndPlotDist(rfit, [], hax(2));
title({'Fisher transformed'; 'corr coeff'})
pbaspect([1.3 1 1])
xlim([-1 1]*0.2)
beautify([12 14 15]);

resizeImageForPub('portrait')

export_fig -r150 -transparent images/tao-monte-carlo-null-hyp-distributions.png

%% test monte carlo error bounds idea with fake time series
% created images/bad-regression.png

xx = linspace(-1,1,5000);
noise = randn(size(xx));

% assume I observe yy0
yy0 = xx + 3*noise;
% do my regression
[coeff, conf, dof, err] = dcregress(xx, yy0, [], 0, 0, 0);
% and get this slope
m1 = coeff(2);
% and intercept
c1 = coeff(1);
[NoiseAmp, NoiseSlope]  = EstimateNoiseSpectrum(yy0-m1*xx-c1, [], 0);

m = nan([5000 1]);
tic;
for mc=1:5000
    % generate noise with same properties
    gennoise = synthetic_timeseries_known_spectrum(length(xx), 1, ...
                                                   NoiseAmp, NoiseSlope)';
    yy = m1*xx + c1 + gennoise;
    [coeff, conf, dof, err] = dcregress(xx, yy, [], 0, 0, 0);
    m(mc) = coeff(2);
    mconf(mc) = conf(2);
end

hold on; histogram(m);
linex(calc95(m));
linex(m1);
linex(1, [], 'r');
linex(mean(m) + [-1 1]*mean(mconf));

%% corr dht, T
[c, lags] = GappyCorrelation(T);
c = abs(c);
mask = lags>0 & c>0;
c = c(mask);
lags = lags(mask);

range = 1:100;
[y0, X] = exp_fit(lags(range), c(range), 1)

Tfill = cut_nan(T);
dhtfill = cut_nan(dht);

dcregress(dhtfill(1:30:end)', Tfill(1:30:end)', [], 0, 1, 0);
dcregress(dht, T, [], 0, 1, 0);

%%

E = [1 1; 1 -1; 1 -2];
y = [1; 2; 4]