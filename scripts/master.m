% get & set options
[opt, plotopt] = DefaultOptions;
plotopt.ploterr = 1;
plotopt.nmode = [1 2];

tao = ReadTaoTriton;
% SaveSlopeBounds(tao, opt);
opt.numMC = 1500;
InferModeShape(opt, tao);
RunTests(opt);

[hax, supax] = PlotModeMap(plotopt);
supax.Title.Visible = 'off';
savepdf('images/02-28-bc2m1.pdf');

[hax, supax] = PlotModeMap(plotopt, 3:11, 2:6);
supax.Title.Visible = 'off';
axes(hax(1)); legend('off');
savepdf('images/02-28-bc2m1-zoom.pdf');

% export_fig -transparent -opengl -r150  -p0.005 -nofontswap -depsc images/01-08-bc2m1.png
