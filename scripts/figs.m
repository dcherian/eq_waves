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
plot(bsxfun(@rdivide, Tmode, [-1 1 -1] .* max(abs(Tmode))), -1*Zmode);
linex(0);
ylim([-2000 0]);
ylabel('Z (m)');
title(['Baroclinic mode shapes at (' num2str(lon) 'E, ' ...
       num2str(lat) 'N)']);
legend('1', '2', '3', 'Location', 'SouthEast');
beautify;
set(gcf, 'Position', [675 223 749 872]);

export_fig images/baroclinic-mode-shapes.png

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
