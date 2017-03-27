function [mbound, rbound, m] = CalcSignificanceBounds(x, y, opt, doplot, hax)
% provided two time series x(t), y(t,z) this function
%   a. filters the two.
%   b. Calculates straight line approximation to their spectrum.
%   c. generate numMC timeseries with the same properties and does the regression.
%   d. returns bounds on regression slope
% This should tell you how large the slope between two variables should be to
% conclude that it is significantly different from 0.

% repurposes TestMC and EstimateNoiseSpectrum
    if ~exist('doplot', 'var'), doplot = 0; end
    if ~exist('hax', 'var'), hax = []; end

    if doplot & isempty(hax), figure; end

    disp('Calculating minimum significant slope');

    if opt.TagainstDHT == 0
        error('CalcSignificanceBounds only works for x=DHT, y=T');
    end

    [xAmp, xSlope] = EstimateNoiseSpectrum(x, []);
    % [yAmp, ySlope] = EstimateNoiseSpectrum(y, []);

    % should have more time instants than z-levels.
    % transpose to make sure y = y(t,z);
    if size(y,1) < size(y,2), y = y'; end

    nz = size(y, 2);
    for zz = 1:nz
        [mdist, m{zz}, r, rdist] = TestMC(xSlope, 1, 0, ...
                                          length(y(:,zz)), y(:,zz), xAmp);

        if doplot
            if isempty(hax)
                subplot(ceil(nz/3), 3, zz);
            else
                axes(hax);
            end
            histogram(m{zz});
            linex(calc95(m{zz}));
            title(num2str(zz));
        end

        mbound(zz) = mean(abs(calc95(m{zz})));
        rbound(zz) = mean(abs(calc95(r)));
    end

    % squeeze faaltu dimension.
    if nz == 1, m = m{1}; end
end