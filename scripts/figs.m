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

fnameh = ['../data/dynht/dyn8s170w_dy.cdf'];
%fnameh = ['../data/temp/t5n180w_dy.cdf'];
dht = double(addnan(squeeze(ncread(fnameh,'DYN_13')),1000))';
%dht = double(addnan(squeeze(ncread(fnameh,'T_20')),100))';

%dht = dht(:,5);
filt.halfdef = 'power';
filt.cutoff = sort(2./[0.135 0.155]);
filt.debugflag = 0;

figure;
hdht = PlotSpectrum(dht);

filt.window = 'rect';
hrect = PlotSpectrum(BandPass(dht, filt));

filt.window = 'gauss';
hgauss = PlotSpectrum(BandPass(dht, filt));

filt.window = 'parzen';
hparzen = PlotSpectrum(BandPass(dht, filt));

filt.window = 'butterworth';
hbutt = PlotSpectrum(BandPass(dht, filt));

%linex(1./filt.cutoff);
legend('raw dyn ht', 'rect filt', 'gauss filt', 'parzen filt', 'butterworth');
title(['5 band smoothed spectrum of dyn ht at (170W, 8S) | [' ...
       num2str(sort(filt.cutoff), '%.2f ') ']']);
beautify;
grid on;
linex(0.15, 'bc2m1', 'k');
linex([0.19 0.11 0.085], [], [1 1 1]*0.5);

export_fig -r300 images/filt-compare-dyn
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

%% plot particular mode


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

load bc2m1.mat

lonidx = [11 8 11];
latidx = [4 4 1];

figure; hold on;
set(gcf, 'Position', [516 121 492 674]);
set(gca, 'XAxisLocation', 'top');

for ii=1:length(lonidx)
    label = sprintf('(%dE, %dN)', modes.lon(lonidx(ii)), modes.lat(latidx(ii)));
    h(ii) = plot(data.Twoa{lonidx(ii), latidx(ii)}, -data.Zwoa, ...
                 'DisplayName', label);
end
ylim([-800 0]);
ylabel('Depth (m)'); xlabel('Temp (C)');
legend('Location', 'SouthEast');
beautify;

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
