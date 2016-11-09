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
legend('bc1', 'bc2', 'bc3', 'Location', 'SouthEast');
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
