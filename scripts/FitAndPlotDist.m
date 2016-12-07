function [] = FitAndPlotDist(in, distname, hax)

    if isprop(in, 'ParameterValues')
        % provided with distribution as input
        fitd = in;
    else
        try
            fitd = fitdist(in', distname);
        catch ME
            fitd = fitdist(in, distname);
        end
    end

    if ~exist('hax', 'var')
        figure;
    else
        axes(hax)
    end

    histogram(fitd.InputData.data, 40, 'Normalization', 'pdf');
    hold on;
    limx = xlim;
    xvec = linspace(limx(1), limx(2), 200);
    plot(xvec, fitd.pdf(xvec), 'k-')
    linex(calc95(fitd.InputData.data))