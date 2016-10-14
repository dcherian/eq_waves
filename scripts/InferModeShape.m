% Get mode shape from data by band pass filtering around expected frequency

% This is how the script works for each point in the TAO array
% 1) Read TAO temp data
% 2) Read dynamic height data and band pass filter around given
%    band
% 3) At each depth, band pass filter interpolated temp data, find missing
%    data in both DHT & temp, and remove them. Then, do the inversion and
%    calculate "inferred mode" (Imode) again at each depth.
% 4) Normalize by *max amplitude* and save data to Imode.mat.

function [modes] = InferModeShape(opt, lonrange, latrange)

  %% Data & Parameters

  ticstart = tic;
  datadir = '../data/';

  modes.lat = [8 5 2 0 -2 -5 -8]; % N
  % modes.lon = [95 110 125 140 137 147 155 156 165 170 180]; % W
  % edited lon values to make neater subplots ...
  % (some lons dont have enough data)
  % values in +ve East.
  modes.lon = -1 * fliplr([95 110 125 140 155 170 180 -165 -156 -147 -137]);

  nlon = length(modes.lon);
  nlat = length(modes.lat);

  if ~exist('lonrange', 'var')
      lonrange = 1:nlon;
  end
  if ~exist('latrange', 'var')
      latrange = 1:nlat;
  end

  %% Iterate!
  clear data
  for mm=lonrange
      for nn=latrange

          if modes.lon(mm) > 0
              lonstr = 'e';
          else
              lonstr = 'w';
          end

          if modes.lat(nn) < 0
              latstr = 's';
          else
              latstr = 'n';
          end

          % TAO filenames
          fnamet = [datadir 'temp/t',   num2str(abs(modes.lat(nn))), ...
                    latstr,num2str(abs(modes.lon(mm))),lonstr,'_dy.cdf'];
          fnameh = [datadir 'dynht/dyn',num2str(abs(modes.lat(nn))), ...
                    latstr,num2str(abs(modes.lon(mm))),lonstr,'_dy.cdf'];

          % make sure observations exist before continuing
          if ~exist(fnamet,'file'), continue; end
          if ~exist(fnameh,'file'), continue; end

          fprintf(['TAO at (%.1f %s, %.1f %s)\n'], ...
                  abs(modes.lon(mm)), upper(lonstr), ...
                  abs(modes.lat(nn)), upper(latstr));

          % read in TAO measurements
          tao.T = double(addnan(squeeze(ncread(fnamet,'T_20')),100));
          modes.depth{mm,nn} = squeeze(ncread(fnamet,'depth'));
          data.depth{mm,nn}  = modes.depth{mm,nn};
          tao.dht = ...
              double(addnan(squeeze(ncread(fnameh,'DYN_13')),1000))';
          tao.dht = tao.dht - nanmean(tao.dht);

          % Read in time vectors for both variables
          tao.timedht = ncread(fnameh, 'time');
          tao.timetemp = ncread(fnamet, 'time');

          % Bandpass filter dynamic height
          dhtfilt = BandPass(tao.dht, opt.filt);

          if opt.debug
              opt.hdbg = figure;
              hdbgax1 = subplot(211);
              PlotSpectrum(tao.dht);
              PlotSpectrum(dhtfilt);
              linex(1./opt.windows);
          end

          % Make sure I'm using same time interval for both variables
          % Needed because dynamic height is not available at all
          % time points with temperature measurement
          start = find(tao.timetemp == tao.timedht(1));
          stop = find(tao.timetemp == tao.timedht(end));
          range = start:stop;
          assert(all(tao.timetemp(range) == tao.timedht));

          % nans
          infer_mode = nan(size(modes.depth{mm,nn}));
          infer_mode_error = nan(size(modes.depth{mm,nn}));

          dof = infer_mode;
          corrcoeff = infer_mode;
          Tstd = infer_mode;

          % iterate over depths and regress at each
          clear Tvec Tfilt Tstd
          for ii = 1:length(modes.depth{mm,nn})
              if opt.filter_temp
                  % band pass filter temperature data
                  Tfilt(ii,:) = BandPass(tao.T(ii,:), opt.filt);
              else
                  Tfilt(ii,:) = tao.T(ii,:);
              end

              if opt.debug
                  axes(hdbgax1);
                  PlotSpectrum(cut_nan(tao.T(ii,:)));
                  PlotSpectrum(cut_nan(Tfilt(ii,:)));
              end

              % save temperature series for reference and calculate std dev
              data.Tfilt{mm,nn,ii} = Tfilt(ii,:);
          end

          Tstd = nanstd(Tfilt(:,range)')';
          [infer_mode, infer_mode_error, corrcoeff, dof] ...
              = DoRegression(dhtfilt, Tfilt(:, range), opt);

          [imnorm, imax] = nanmax(abs(infer_mode(:,1)));
          modes.InferredModeOLS{mm,nn} = infer_mode(:,1)./imnorm;
          modes.InferredModeErrorOLS{mm,nn} = ...
              infer_mode_error(:,1)./imnorm;

          [imnorm, imax] = nanmax(abs(infer_mode(:,2)));
          modes.InferredModeWTLS{mm,nn} = infer_mode(:,2)./imnorm;
          modes.InferredModeErrorWTLS{mm,nn} = ...
              infer_mode_error(:,2)./imnorm;

          modes.dof{mm,nn} = dof;
          modes.corr{mm,nn} = corrcoeff;

          modes.timedht{mm,nn} = tao.timedht;
          modes.dhtfilt{mm,nn} = dhtfilt;
          modes.Tstd{mm,nn} = Tstd;
      end
  end

  hash = githash([mfilename('fullpath') '.m']);

  modes.comment = ['InferredMode(modes.lon,modes.lat,modes.depth)' ...
                   'is the ' ...
                   'mode structure inferred by regressing ' ...
                   'band passed dyn. ht against band passed' ...
                   'temperature. ' ...
                  'Tfilt = bandpassed temperature (cell array) | ' ...
                  'dhtfilt = band passed dynamic height | ' ...
                  'Tstd = standard dev of temp time series at depth'];

  % don't overwrite if inferring at all locations
  if isequal(lonrange, 1:nlon) & isequal(latrange, 1:nlat)
      save([opt.name '-' opt.filt.window '.mat'], ...
           'modes', 'data', 'opt', 'hash');
  end
  toc(ticstart);
end
