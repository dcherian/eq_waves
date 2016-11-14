%% typical mode shapes

load('./flat-bot-modes.mat');

mm = 11;
nn = 4;

linestylemode = {'--', '-', '-.'}; % line style for theoretical modes.

figure;
clf; hold on;
for ii=1:3
    mode = squeeze(flatbot.IdealTempMode(mm, nn, :, ii));
    plot(mode./max(mode), -1 * flatbot.zTmode, ...
         'Color', 'k', 'LineStyle', linestylemode{ii})
end
hax = gca;
hax.XAxisLocation = 'top';
ylim([-1000 0])
xlim([-0.25 1]);
yticks = hax.YTick;
liney([-1 -25 -50 -75 -100 -125 -150 -200 -250 -300 -500 -750]);
hax.YTick = yticks;
legend('T_{bc1}', 'T_{bc2}', 'T_{bc3}', 'Location', 'SouthEast');
linex(0);
resizeImageForPub('onecolumn');
hax.Position(2) = 0.05;
pbaspect([1 1.615 1]);
beautify([13 14 15], 'Times');
xlabel({'Flat-bottom temp mode'; ...
        ['amplitudes at ' getTitleString(flatbot.lon(mm), ...
                                         flatbot.lat(nn))]});
ylabel('Z (m)');

export_fig -c[Inf,0,0,0] -despc2 images/flat-bottom-modes.pdf

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
