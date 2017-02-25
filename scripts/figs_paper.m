%% more mode shapes
load('./flat-bot-modes.mat');
label = 'abcdefg';

hfig = figure;
for ii=1:8
    hax(ii) = subplot(2,4,ii);
end
hax(4) = subplot(2,4, [4 8]);
hax(6).delete;
resizeImageForPub('portrait');
hfig.Position(4) = 600;

% typical mode shapes
mm = 11;
nn = 4;

linestylemode = {'--', '-', '-.'}; % line style for theoretical modes.
flatbot.IdealWMode(mm,nn,:,1) = flatbot.IdealWMode(mm,nn,:,1) * -1;

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

for xx = 1:3
    axes(hax(xx))
    hax(xx).Position(3) = 0.16;
    hax(xx).XAxisLocation = 'top';
    ylim([-1000 0])
    xlim([-0.25 1]);
    yticks = hax(xx).YTick;
    hly = liney([-1 -25 -50 -75 -100 -125 -150 -200 -250 -300 -500 -750]);
    hax(xx).YTick = yticks;
    hx0 = linex(0);
    uistack([hly hx0], 'bottom');
    pbaspect([1 2.5 1]);
    if xx == 1
        hleg2 = legend('$\bar{T}$', '$d\bar{T}/dz$',...
                       'Location', 'SouthEast');
        xlabel(['(a) Mean state']);
        hleg2.Position = [0.230 0.5817 0.1333 0.0575];
    else
        hax(xx).YTickLabels = {};
    end
    if xx == 2
        legend('off');
        %hleg(2) = legend('w_{bc1}', 'w_{bc2}', 'w_{bc3}', 'Location', ...
        %'SouthEast');
        xlabel(['(b) $w$ mode']);
    end
    if xx == 3
        hleg4 = legend('bc1', 'bc2', 'bc3', 'Location', ...
                       'SouthEast');
        xlabel(['(c) $T$ mode']);
        hleg4.Position = [0.4 0.5800 0.1111 0.0833];
    end
    beautify([13 14 15], 'Times');
end

hax(1).YLabel.String = 'z (m)';
hleg2.Interpreter = 'latex';
hax(2).XLabel.Interpreter = 'latex';
hax(3).XLabel.Interpreter = 'latex';

% inferred mode with phase lag plot
mm = 5; nn = 4;
[~, plotopt] = DefaultOptions;
plotopt.plotPhaseLag = 1;
plotopt.MarkWaterDepth = 0;

PlotMode('bc2m1-butterworth', mm, nn, plotopt, hax(5));
legend('off');
xlabel(['(d) ' hax(5).Title.String])
hax(5).Title.delete;
beautify([12 13 14], 'Times');

mm = 10; nn = 4;
PlotMode('bc2m1-butterworth', mm, nn, plotopt, hax(7));
xlabel(['(e) ' hax(7).Title.String])
hax(7).Title.delete;
hleg5 = legend;
hleg5.Position = [0.2482    0.100    0.2389    0.1258];
hleg5.String{1} = ['T_{infer}'];
beautify([12 13 14], 'Times');

linkaxes(hax([5 7]), 'xy')
hax(5).YLim = [-600 0];
hax(5).XLim = [-0.5 1.2]*1;
hax(5).YLabel.String = 'z (m)';

hax(5).Position(2) = 0.07;
hax(7).Position(2) = 0.07;
hax(5).Position(3) = 0.24;
hax(7).Position(3) = 0.24;
hax(7).Position(1) = 0.47;
hax(7).YTickLabels = {};

htxt.delete;
htxt = text(-0.66, 1.25, 'Inferred mode shapes', ...
            'FontSize', 16, 'FontName', 'Times', ...
            'Units', 'normalized');

% deep flat-bottom temperature mode
axes(hax(4))
hax(4).Position(1) = 0.80;
hold on;
n_mode = 2;
ilon = [5 11];
ilat = [4 4];
colors = get(groot, 'DefaultAxesColorOrder');
for ii=1:length(ilon)
    mm = ilon(ii);
    nn = ilat(ii);
    mode = squeeze(flatbot.IdealTempMode(mm, nn, :, n_mode));
    plot(mode./max(mode), -1 * flatbot.zTmode, ...
         'DisplayName', ...
         ['T-mode ' getTitleString(flatbot.lon(mm),  flatbot.lat(nn))], ...
         'Color', colors(ii,:));
    mode = squeeze(flatbot.IdealUMode(mm, nn, :, n_mode));
    plot(mode./max(abs(mode)), -1 * flatbot.zTmode, ...
         'DisplayName', ...
         ['u-mode ' getTitleString(flatbot.lon(mm), flatbot.lat(nn))], ...
         'LineStyle', '--', 'Color', colors(ii,:));
end
hleg1 = legend('Location', 'SouthEast');
linex(0);
ylim([-3500 0]);
xlim([-1 1]);
beautify([12 13 14], 'Times');
title('(f) T_{bc2}')
ylabel('z (m)')
hax(4).XAxisLocation = 'top';
pbaspect([1 4 1]);
hleg1.Position(2) = 0.10;

savepdf('images/multi-panel-modes.pdf');
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