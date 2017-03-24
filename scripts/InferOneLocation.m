function [] = InferOneLocation(mm, nn, opt, plotopt)
   if ~exist('opt', 'var') | isempty(opt)
       [opt, ~] = DefaultOptions;
   end
   if ~exist('plotopt', 'var') | isempty(plotopt)
       [~, plotopt] = DefaultOptions;
   end

   tao = ReadTaoTriton(mm,nn);
   modes = InferModeShape(opt, tao, mm, nn);

   zz = 5;
   figure;
   hax(1) = subplot(121);
   PlotMode(modes, mm, nn, plotopt, hax(1));

   hax(2) = subplot(122);
   keyboard;
   % PlotSpectrum(EstimateNoiseSpectrum(tao.T{mm,nn}(zz,:), opt, 1, hax(2)))
   PlotSpectrum(tao.dht{mm,nn});
   PlotSpectrum(BandPass(tao.dht{mm,nn}, opt.filt));
   % PlotSpectrum(tao.T{mm,nn}(zz,:))
   % PlotSpectrum(BandPass(tao.T{mm,nn}(zz,:), opt.filt));
end