
%% figure out data frequency slopes
tao = ReadTaoTriton(11,4);

figure;
subplot(121)
hspec = PlotSpectrum(tao.dht{11,4});
PlotSpectrum(synthetic_timeseries_known_spectrum(length(tao.dht{11,4}),1,0.05,-2.4));
uistack(hspec, 'top');

subplot(122);
hspec = PlotSpectrum(tao.dht{11,4});
PlotSpectrum(synthetic_timeseries_known_spectrum(length(tao.dht{11,4}),1,0.05,-2));
uistack(hspec, 'top');

%% Monte Carlo regression Slope

%slope = 0.15;
%intercept = 0.2;

spectralSlope = -2;
len = 2000; % number of data points

xx = synthetic_timeseries_known_spectrum(len, 1, 1, spectralSlope);
opt = DefaultOptions;
xx = BandPass(xx, opt.filt);
xx = cut_nan(xx - nanmean(xx));

tic;
numMC = 1e3;
clear c m r
for ii=1:numMC
    % xnoise = 0;
    % yy = slope * (xx + xnoise) + intercept;
    % ynoise = 0.2*max(abs(yy))*whitenoise(size(yy));
    % yy = yy + ynoise;

    yy = synthetic_timeseries_known_spectrum(len, ...
                                             1, 2, spectralSlope);
    yy = BandPass(yy, opt.filt);
    yy = cut_nan(yy - nanmean(yy));

    coeff = dcregress(xx, yy, [], 0, 0, 0);
    c(ii) = coeff(1);
    m(ii) = coeff(2);

    rmat = corrcoef(cut_nan(xx), cut_nan(yy));
    r(ii)= rmat(1,2);
end
toc;

%[coeff, conf, dof] = dcregress(xx, yy, length(yy)-2, 0, 0);
% disp(['actual: slope = ' num2str(slope, '%.3f') ' | intercept = ' ...
%       num2str(intercept, '%.3f')]);
% disp(['dcregress: slope = [' num2str([coeff(2)-conf(2) coeff(2)+conf(2)], '%.3f ') ...
%       '] | intercept = [' num2str([coeff(1)-conf(1) coeff(1)+conf(1)], ...
%                                   '%.3f ') ']']);

nbins = 40;
figure;
subplot(221)
plot(xx,yy, '*');
title(['Spectral slope = ' num2str(spectralSlope)]);
beautify;

subplot(222)
histogram(r, nbins, 'Normalization', 'countdensity', ...
          'EdgeColor', 'none');
linex(calc95(r), 'r');
title('corr. coeff');
beautify;

subplot(223);
hsl = histogram(m, nbins, 'Normalization', 'countdensity', ...
                'EdgeColor', 'none');
xbin = hsl.BinEdges;
%linex(slope, [], 'k')
%linex([coeff(2) coeff(2)-conf(2) coeff(2)+conf(2)]);
linex(calc95(m), [], 'r');
title('slope');
beautify;

subplot(224);
histogram(c, nbins, 'Normalization', 'countdensity', 'EdgeColor', 'none');
%linex(intercept, [], 'k');
%linex([coeff(1) coeff(1)-conf(1) coeff(1)+conf(1)]);
linex(calc95(c), [], 'r');
title('intercept');
beautify;

resizeImageForPub('portrait');
export_fig images/monte-carlo-regression--2.png
