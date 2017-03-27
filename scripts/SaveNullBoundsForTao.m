function [] = SaveNullBoundsForTao(tao, opt, lmm, lnn)

    if ~exist('lmm', 'var') | isempty(lmm), lmm = 1:11; end
    if ~exist('lnn', 'var') | isempty(lnn), lnn = 1:7; end

    if exist('bounds.mat', 'file')
        load bounds.mat
    else
        mbound = cell(size(tao.dht, 1), size(tao.dht, 2));
        rbound = cell(size(tao.dht, 1), size(tao.dht, 2));
        mdist = cell(size(tao.dht, 1), size(tao.dht, 2));
    end

    hash = githash([mfilename('fullpath') '.m']);

    tic;
    disp('Monte carlo estimation of null bounds on regression slope');
    for mm=lmm
        for nn=lnn
            if isempty(tao.dht{mm,nn})
                mbound{mm,nn} = NaN;
                rbound{mm,nn} = NaN;
                continue;
            end

            disp([mm nn]);
            for zz=1:length(tao.depth{mm,nn})
                if isempty(tao.T{mm,nn}(zz,:))
                    mbound{mm,nn}(zz) = NaN;
                    rbound{mm,nn}(zz) = NaN;
                    continue;
                end

                range = findCommonTimeRange(tao.timedht{mm,nn}, ...
                                            tao.timetemp{mm,nn});
                [mbound{mm,nn}(zz), rbound{mm,nn}(zz), mdist{mm,nn}{zz}] = ...
                    CalcSignificanceBounds(tao.dht{mm,nn}, ...
                                           tao.T{mm,nn}(zz,range), ...
                                           opt);
            end

            save('./bounds.mat', 'mbound', 'rbound', 'mdist', 'hash');
        end
    end
    toc;
end