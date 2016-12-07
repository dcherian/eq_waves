% get & set options
[opt, plotopt] = DefaultOptions;
plotopt.plotWTLS = 0;
plotopt.ploterr = 1;
plotopt.nmode = [1 2];

tao = ReadTaoTriton;
InferModeShape(opt, tao);
PlotModeMap(plotopt);

export_fig  -p0.005 -nofontswap -depsc images/12-07-bc2m1.pdf
