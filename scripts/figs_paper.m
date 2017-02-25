%% more mode shapes
load('./flat-bot-modes.mat');

figure;
hold on;
n_mode = 2;
ilon = [5 11];
ilat = [4 4];
for ii=1:length(ilon)
    mm = ilon(ii);
    nn = ilat(ii);
    mode = squeeze(flatbot.IdealTempMode(mm, nn, :, n_mode));
    plot(mode./max(mode), -1 * flatbot.zTmode, ...
         'DisplayName', getTitleString(flatbot.lon(mm), flatbot.lat(nn)));
end
legend('Location', 'SouthEast');
linex(0);
ylim([-2500 0]);
xlim([-0.15 1]);
beautify([12 13 14], 'Times');
resizeImageForPub('onecolumn');
title('Mode-2 temperature mode shapes')
ylabel('Depth (m)')
hax = gca;
hax.XAxisLocation = 'top';

export_fig -c[Inf,0,Inf,0] images/flat-bottom-modes-deeper.pdf

%% typical mode shapes

load('./flat-bot-modes.mat');

mm = 11;
nn = 4;

linestylemode = {'--', '-', '-.'}; % line style for theoretical modes.

flatbot.IdealWMode(mm,nn,:,1) = flatbot.IdealWMode(mm,nn,:,1) * -1;

figure;
hax = packfig(1,3);
hax(1).NextPlot = 'add';
hax(2).NextPlot = 'add';
hax(3).NextPlot = 'add';

axes(hax(1));
plot(flatbot.Twoa{mm,nn}./max(flatbot.Twoa{mm,nn}), -1 * flatbot.Zwoa);
plot(hax(1), squeeze(flatbot.dTdz(mm,nn,:)) ...
     ./max(squeeze(flatbot.dTdz(mm,nn,:))), -1*flatbot.zTmode, 'k-');

for ii=1:3
    axes(hax(2));
    mode = squeeze(flatbot.IdealWMode(mm, nn, :, ii));
    plot(-1 * mode./max(abs(mode)), -1 * flatbot.zTmode, ...
         'Color', 'k', 'LineStyle', linestylemode{ii});

    axes(hax(3));
    mode = squeeze(flatbot.IdealTempMode(mm, nn, :, ii));
    plot(mode./max(mode), -1 * flatbot.zTmode, ...
         'Color', 'k', 'LineStyle', linestylemode{ii})
end

label = 'abc';
for xx = 1:length(hax)
    axes(hax(xx))
    title(['(' label(xx) ')']);
    hax(xx).XAxisLocation = 'top';
    ylim([-1000 0])
    xlim([-0.25 1]);
    yticks = hax(xx).YTick;
    hly = liney([-1 -25 -50 -75 -100 -125 -150 -200 -250 -300 -500 -750]);
    hax(xx).YTick = yticks;
    hx0 = linex(0);
    uistack([hly hx0], 'bottom');
    pbaspect([1 1.615 1]);
    if xx == 1
        hleg(1) = legend('T', 'dT/dz', 'Location', 'SouthEast');
    end
    if xx == 2
        hleg(2) = legend('w_{bc1}', 'w_{bc2}', 'w_{bc3}', 'Location', ...
                         'SouthEast');
    end
    if xx == 3
        hleg(3) = legend('T_{bc1}', 'T_{bc2}', 'T_{bc3}', 'Location', ...
                         'SouthEast');
    end
    beautify([13 14 15], 'Times');
end

resizeImageForPub('portrait');
hax(1).YLabel.String = 'Z (m)';
[supax] = suplabel(['Flat-bottom w and T mode amplitudes at ' ...
                    getTitleString(flatbot.lon(mm), flatbot.lat(nn))], ...
                     't');
axes(supax); beautify([13 14 15], 'Times');
supax.Position(4) = 0.8;
hleg(3).Position(1) = 0.8;

export_fig -c[Inf,0,Inf,0] -despc2 images/flat-bottom-modes.pdf

%% dynht spectrum with filter bounds

[opt,~] = DefaultOptions;
hf = openfig('../images/dyn_ht_5S_5N.fig');
hax = hf.Children;
% hax(1) is main spectrum
% hax(2) is colorbar

axes(hax(1));
handles = liney(2./opt.filt.cutoff);
title({'Farrar & Durland (2012) spectrum'; '(5S-5N)'});
beautify([13 14 15], 'Times');
caxis([1 2.5]);
colormap(cbrewer('seq', 'OrRd', 20));
hcb = colorbar;
hcb.Label.String = 'Log_{10} of spectral density (cm^2/cpd/deg^{-1})';
hcb.Label.FontSize = 13;

hax(2).Color = 'None';
hax(2).YAxis.Color = 'None';
hax(2).XAxis.Color = 'None';
hax(2).Children(1).FontSize = 12;
hax(2).Children(1).FontName = 'Times';
hax(2).Children(1).Color = hcb.Label.Color;
hax(2).Children(end).delete;
for ii=2:5
    hax(2).Children(ii).LineWidth = 1;
    hax(2).Children(ii).Color = hcb.Label.Color;
end

resizeImageForPub('onecolumn');

hax(1).Position(1) = 0.19;
hax(2).Position(1) = 0.88;
hax(2).Position(4) = 0.7;

hax(1).Color = 'None';

export_fig -c[Inf,0,Inf,0] -despc ./images/farrar-durland-spectrum.pdf

%% background noise estimation

opt = DefaultOptions;

hfig = figure;
hax(1) = subplot(211);
hax(2) = subplot(212);
kk = 1;

ilon = [4 11];
ilat = [7 7];
for ii = 1:length(ilon)
    mm = ilon(ii);
    nn = ilat(ii);
    data = ReadTaoTriton(mm,nn);
    EstimateNoiseSpectrum(data.dht{mm,nn}, opt, 1, hax(kk));
    hax(kk).Title.String = getTitleString(data.lon(mm),data.lat(nn));
    %EstimateNoiseSpectrum(data.dht{mm,nn}, opt, 1, hax(kk+2));
    if kk == 1
        legend('off');
        xlabel('');
    else
        hleg = legend;
        hleg.String{1} = 'Dyn. ht.';
        hleg.String{2} = 'Filtered dyn. ht.';
        hleg.Position = [ 0.5741    0.2905    0.4407    0.1810];
    end
        axes(hax(kk)); beautify([12 13 14], 'Times');
        kk = kk+1;
end
linkaxes(hax, 'xy');
hax(1).XLim = [0.05 0.5];
hax(1).YLim = [1e-2 1e4];
resizeImageForPub('onecolumn');

% set(hfig, 'Position', [300 415 540 285]);
hfig.Position = [130 415 270 420];
% hleg.Position =  [0.1256 0.3040 0.2286 0.2667];
savepdf('images/estimate-noise-spectrum.pdf')

%% bathy map

datadir = '../data/';
if ~exist('etopo', 'var')
    etoponame = ['../data/ETOPO2v2g_f4.nc4'];
    disp(' Loading ETOPO2v2 data.');
    etopo.x = ncread(etoponame, 'x');
    etopo.y = ncread(etoponame, 'y');
    etopo.z = ncread(etoponame, 'z');
end

tao.lat = [8 5 2 0 -2 -5 -8];
% values in +ve East.
tao.lon = -1 * fliplr([95 110 125 140 155 170 180 ...
                    -165 -156 -147 -137]);

figure; hold on;
PlotBathy(etopo, [-180 -80], [-1 1]*10, 10);
PlotBathy(etopo, [+110-360 -180], [-1 1]*10, 10);
for mm=1:11
    for nn=1:7
        if tao.lon(mm) > 0
            lonstr = 'e';
        else
            lonstr = 'w';
        end

        if tao.lat(nn) < 0
            latstr = 's';
        else
            latstr = 'n';
        end

        % TAO filenames
        fnamet = [datadir 'temp/t',   num2str(abs(tao.lat(nn))), ...
                  latstr,num2str(abs(tao.lon(mm))),lonstr,'_dy.cdf'];
        fnameh = [datadir 'dynht/dyn',num2str(abs(tao.lat(nn))), ...
                  latstr,num2str(abs(tao.lon(mm))),lonstr,'_dy.cdf'];

        % make sure observations exist before continuing
        if ~exist(fnamet,'file'), continue; end
        if ~exist(fnameh,'file'), continue; end

        if tao.lon(mm) > 0
            lon = -360+tao.lon(mm);
        else
            lon = tao.lon(mm);
        end

        h(mm,nn) = plot(lon, tao.lat(nn), 'r.', ...
                        'MarkerSize', 20);
    end
end

heq = liney(0);
uistack(heq, 'bottom');
resizeImageForPub('portrait')
pbaspect([4 1 1]);
xlabel('Longitude');
ylabel('Latitude');
beautify([12 12 14], 'Times'); box on;

hax = gca;
ticks = hax.XTick;
for ii=1:length(hax.XTick)
    if ticks(ii) < -180
        hax.XTickLabel{ii} = [num2str(360+ticks(ii)) '^oE'];
    else
        hax.XTickLabel{ii} = [num2str(-ticks(ii)) '^oW'];
    end
end
hax.YAxis.TickLabelFormat = '%g^o';

export_fig -c[Inf,0,Inf,0] -despc2 images/bathy-mooring-locations.pdf

%% phase lag plot

mm = 5; nn = 4;
[~, plotopt] = DefaultOptions;
plotopt.plotPhaseLag = 1;
plotopt.MarkWaterDepth = 0;

figure;
hax = packfig(1,2);
PlotMode('bc2m1-butterworth', mm, nn, plotopt, hax(1));
legend('off');
beautify([12 13 14], 'Times');

mm = 10; nn = 4;
PlotMode('bc2m1-butterworth', mm, nn, plotopt, hax(2));
beautify([12 13 14], 'Times');

linkaxes(hax, 'xy');
hax(1).YLim = [-600 0];
hax(1).XLim = [-0.5 1.2]*1;
hax(1).YLabel.String = 'Depth (m)';

pbaspect(hax(1), [1 1.412 1]);
pbaspect(hax(2), [1 1.412 1]);

set(gcf, 'Position', [360   325   540   373]);
hax(1).Position(2) = 0.06;
hax(2).Position(2) = 0.06;
axes(hax(2)); hleg = legend; hleg.Position(1) = 0.7;
resizeImageForPub('portrait')
savepdf('images/phase-lag.pdf');