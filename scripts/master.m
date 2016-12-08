% get & set options
[opt, plotopt] = DefaultOptions;
plotopt.plotWTLS = 0;
plotopt.ploterr = 1;
plotopt.nmode = [1 2];

tao = ReadTaoTriton;
InferModeShape(opt, tao);
RunTests(opt);
[hax, supax] = PlotModeMap(plotopt);
supax.Title.Visible = 'off';

export_fig -transparent -opengl -r150  -p0.005 -nofontswap -depsc images/12-08-bc2m1.png
