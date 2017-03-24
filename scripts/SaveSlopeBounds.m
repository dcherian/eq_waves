function [] = SaveNullBounds(tao, opt)

    tic;
    disp('Monte carlo estimation of null bounds on regression slope');
    for mm=1:11
        for nn=1:7
            if isempty(tao.dht{mm,nn})
                mbound{mm,nn} = NaN;
                rbound{mm,nn} = NaN;
            end

            for zz=1:length(tao.depth{mm,nn})
                if isempty(tao.T{mm,nn}(zz,:))
                    mbound{mm,nn}(zz) = NaN;
                    rbound{mm,nn}(zz) = NaN;
                    continue;
                end

                range = findCommonTimeRange(tao.timedht{mm,nn}, ...
                                            tao.timetemp{mm,nn});
                [mbound{mm,nn}(zz), rbound{mm,nn}(zz)] = ...
                    CalcSignificanceBounds(tao.dht{mm,nn}, ...
                                           tao.T{mm,nn}(zz,range));
            end
        end
    end
    toc;
    githash([mfilename('fullpath') '.m']);

    save('./bounds.mat', 'mbound', 'rbound');
    toc;
end

function [mbound, rbound] = CalcSignificanceBounds(x, y, opt, doplot)
% provided two time series x,y this function
%   a. filters the two.
%   b. Calculates straight line approximation to their spectrum.
%   c. generate numMC timeseries with the same properties and does the regression.
%   d. returns bounds on regression slope
% This should tell you how large the slope between two variables should be to
% conclude that it is significantly different from 0.

% repurposes TestMC and EstimateNoiseSpectrum
    if ~exist('doplot', 'var'), doplot = 0; end

    numMC = 1000;
    [yAmp, ySlope] = EstimateNoiseSpectrum(y, []);

    [mdist, m, r, rdist] = TestMC(ySlope, 1, 0, length(y), x, yAmp);

    if doplot
        figure;
        histogram(m);
        linex(calc95(m));
    end

    mbound = mean(abs(calc95(m)));
    rbound = mean(abs(calc95(r)));
end