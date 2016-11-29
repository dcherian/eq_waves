% Reads TAO/TRITON data and returns structure tao.
function [tao] = ReadTaoTriton(lonrange, latrange)

    ticstart = tic;
    disp('Reading TAO/TRITON data');

    datadir = '../data/';

    tao.name = 'tao';
    tao.lat = [8 5 2 0 -2 -5 -8];
    % modes.lon = [95 110 125 140 137 147 155 156 165 170 180]; % W
    % edited lon values to make neater subplots ...
    % (some lons dont have enough data)
    % values in +ve East.
    tao.lon = -1 * fliplr([95 110 125 140 155 170 180 ...
                        -165 -156 -147 -137]);

    nlon = length(tao.lon);
    nlat = length(tao.lat);

    if ~exist('lonrange', 'var')
        lonrange = 1:nlon;
    end
    if ~exist('latrange', 'var')
        latrange = 1:nlat;
    end

    for mm=lonrange
        for nn=latrange
            if tao.lon(mm) > 0
                lonstr = 'e';
            else
                lonstr = 'w';
            end

            if tao.lat(nn) < 0
                latstr = 's';
            else
                latstr = 'n';
            end

            % TAO filenames
            fnamet = [datadir 'temp/t',   num2str(abs(tao.lat(nn))), ...
                      latstr,num2str(abs(tao.lon(mm))),lonstr,'_dy.cdf'];
            fnameh = [datadir 'dynht/dyn',num2str(abs(tao.lat(nn))), ...
                      latstr,num2str(abs(tao.lon(mm))),lonstr,'_dy.cdf'];

            % make sure observations exist before continuing
            if ~exist(fnamet,'file'), continue; end
            if ~exist(fnameh,'file'), continue; end

            fprintf(['TAO/TRITON at (%.1f %s, %.1f %s)\n'], ...
                    abs(tao.lon(mm)), upper(lonstr), ...
                    abs(tao.lat(nn)), upper(latstr));

            % read in TAO measurements
            tao.T{mm,nn} = double(addnan(squeeze(ncread(fnamet,'T_20')),100));
            tao.depth{mm,nn} = squeeze(ncread(fnamet,'depth'));
            tao.dht{mm,nn} = ...
                double(addnan(squeeze(ncread(fnameh,'DYN_13')),1000))';

            % Read in time vectors for both variables
            tao.timedht{mm,nn} = ncread(fnameh, 'time');
            tao.timetemp{mm,nn} = ncread(fnamet, 'time');
        end
    end

    toc(ticstart);
end