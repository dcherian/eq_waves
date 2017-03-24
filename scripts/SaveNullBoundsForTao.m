function [] = SaveNullBoundsForTao(tao, opt)

    hash = githash([mfilename('fullpath') '.m']);

    mbound = cell(size(tao.dht, 1), size(tao.dht, 2));
    rbound = cell(size(tao.dht, 1), size(tao.dht, 2));

    tic;
    disp('Monte carlo estimation of null bounds on regression slope');
    for mm=1:11
        for nn=1:7
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
                [mbound{mm,nn}(zz), rbound{mm,nn}(zz)] = ...
                    CalcSignificanceBounds(tao.dht{mm,nn}, ...
                                           tao.T{mm,nn}(zz,range));
            end

            save('./bounds.mat', 'mbound', 'rbound', 'hash');
        end
    end
    toc;
end