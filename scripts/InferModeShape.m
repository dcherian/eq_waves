% Get mode shape from data by band pass filtering around expected frequency

% This is how the script works for each point in the TAO array
% 1) Load WOA data and calculate theoretical 0 mean vel. + flat
%    bottom mode shapes. The T mode is obtained by multiplying the
%    vertical mode by dT/dz (from WOA).
% 2) Read TAO temp data
% 3) Read dynamic height data and band pass filter around given
%    band
% 4) At each depth, band pass filter interpolated temp data, find missing
%    data in both DHT & temp, and remove them. Then, do the inversion and
%    calculate "inferred mode" (Imode) again at each depth.
% 5) Normalize by *max amplitude* and save data to Imode.mat.

function [] = InferModeShape(opt)

  %% Data & Parameters

  datadir = '../data/';

  if exist('temp','var') ~= 1 && exist('sal','var') ~=1
      fprintf('\n Loading WOA 05 data \n\n');
      woa = load([datadir 'woa05.mat']);
  end

  modes.lat = [8 5 2 0 -2 -5 -8]; % N
  % modes.lon = [95 110 125 140 137 147 155 156 165 170 180]; % W
  % edited lon values to make neater subplots ...
  % (some lons dont have enough data)
  modes.lon = fliplr([95 110 125 140 155 170 180]);

  nlon = length(modes.lon);
  nlat = length(modes.lat);

  %% Iterate!
  clear data
  for mm=1:nlon
      for nn=1:nlat
          %% Step 1. Calculate theoretical modes

          % Locate (modes.lon,lat)
          ilon = find_approx(woa.X,360-modes.lon(mm),1);
          ilat = find_approx(woa.Y,modes.lat(nn),1);

          % Verify location
          fprintf('\n Actual location = (%.1f W,%.1f N) \n', ...
                  360-woa.X(ilon), woa.Y(ilat));

          % Interpolate to denser grid.
          Zmode = avg1(woa.Z); % Z from woa05.mat
          T = woa.temp(:,ilat,ilon);
          S = woa.sal(:,ilat,ilon);
          dtdz = diff(T)./diff(woa.Z);
          %P = sw_pres(modes.depth,modes.lat);

          N2 = bfrq(S,T,woa.Z,modes.lat(nn));%10^(-6)*ones(32,1);
          [Vmode, Hmode, c] = vertmode(N2,woa.Z,opt.n_modes,0);
          % calculate temperature mode shape
          Tmode = Vmode .* repmat(dtdz,1,size(Vmode,2));

          %% Filter out TAO data

          % clear tbuoy depth dht Tfilt Tfilt1 dhtfilt

          % read in TAO measurements
          if modes.lat(nn) < 0
              fnamet = [datadir 'temp/t',   num2str(abs(modes.lat(nn))), ...
                        's',num2str(modes.lon(mm)),'w_dy.cdf'];
              fnameh = [datadir 'dynht/dyn',num2str(abs(modes.lat(nn))), ...
                        's',num2str(modes.lon(mm)),'w_dy.cdf'];
          else
              fnamet = [datadir 'temp/t',   num2str(abs(modes.lat(nn))), ...
                        'n',num2str(modes.lon(mm)),'w_dy.cdf'];
              fnameh = [datadir 'dynht/dyn',num2str(abs(modes.lat(nn))), ...
                        'n',num2str(modes.lon(mm)),'w_dy.cdf'];
          end

          % Read & interpolate temp
          if ~exist(fnamet,'file'), continue; end
          tao.T = double(addnan(squeeze(ncread(fnamet,'T_20')),100));
          modes.depth{mm,nn} = squeeze(ncread(fnamet,'depth'));
          data.depth{mm,nn}  = modes.depth{mm,nn};

          % interpolate buoy profiles onto standard depth
          % why do this?
          % tao.T = interp1(depth,tao.T,modes.depth);

          % Read dynamic ht
          if ~exist(fnameh,'file'), continue; end
          tao.dht = ...
              double(addnan(squeeze(ncread(fnameh,'DYN_13')),1000))';
          tao.dht = tao.dht - nanmean(tao.dht);

          % Read in time vectors for both variables
          tao.timedht = ncread(fnameh, 'time');
          tao.timetemp = ncread(fnamet, 'time');

          % Bandpass filter dynamic height
          if opt.butterworth
              [b,a] = butter(12, sort(2*pi./opt.windows/(2*pi/2)), 'bandpass');
              dhtfilt = filter(b,a,tao.dht);
          else
              % Running average windows(1) , windows(2) day and subtract
              dhtfilt = conv_band_pass(tao.dht,opt.windows);
          end

          if opt.debug
              hdbg = figure;
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
          Tstd = infer_mode;

          % iterate over depths and regress at each
          clear Tvec Tfilt Tstd
          for ii = 1:length(modes.depth{mm,nn})
              if opt.filter_temp
                  % band pass filter temperature data
                  if opt.butterworth
                      % filter removing NaNs
                      Tvec = cut_nan(tao.T(ii,:));
                      Tfilt(ii,~isnan(tao.T(ii,:))) = ...
                          filter(b,a,Tvec-mean(Tvec));
                      Tfilt(ii,isnan(tao.T(ii,:))) = NaN;
                      filtrange = range;
                  else
                      Tfilt(ii,:) = conv_band_pass(tao.T(ii,range), ...
                                                   opt.windows);
                      filtrange = 1:size(Tfilt,2);
                  end
              else
                  Tfilt(ii,:) = tao.T(ii,:);
              end

              if opt.debug
                  axes(hdbgax1);
                  PlotSpectrum(cut_nan(tao.T(ii,:)));
                  PlotSpectrum(cut_nan(Tfilt(ii,:)));
              end

              % save temperature series for reference and calculate std dev
              data.Tfilt{mm,nn,ii,:} = Tfilt(ii,:);
              Tstd(ii) = nanstd(Tfilt(ii,:));

              % find all nan's in both datasets
              Treduced = Tfilt(ii, filtrange);
              mask = ~(isnan(dhtfilt) | isnan(Treduced));

              % regress to find mode shape
              infer_mode(ii) = dhtfilt(mask)' \ Treduced(mask)';

              if opt.debug
                  figure(hdbg);
                  hdbgax2 = subplot(212);
                  plot(dhtfilt(mask)); hold on;
                  plot(Treduced(mask));
                  title('filtered time series for regression');
                  keyboard;
              end

              %Imode = fill_gap(dhtfilt(filtrange)','linear',15)\fill_gap(Tfilt(:,filtrange)','linear',15);

              modes.InferredMode{mm,nn} = infer_mode./max(abs(infer_mode));
              modes.IdealTempMode(mm,nn,:) = Tmode(:,opt.n_mode);
              data.Tstd{mm,nn} = Tstd;
          end
      end
  end

  modes.zTmode = Zmode;

  hash = githash([mfilename('fullpath') '.m']);
  data.dhtfilt = dhtfilt;
  data.comment = ['Tfilt = bandpassed temperature (cell array) | ' ...
                  'dhtfilt = band passed dynamic height' ...
                  'Tstd = standard dev of temp time series at depth'];

  modes.comment = ['InferredMode(modes.lon,modes.lat,modes.depth)' ...
                   'is the ' ...
                   'mode structure inferred by regressing ' ...
                   'band passed dyn. ht against band passed' ...
                   'temperature. ' ...
                   'modes.IdealTempMode is the theoretical' ...
                   ' temperature mode on grid modes.zTmode'];

  save([opt.name '.mat'], 'modes', 'data', 'opt', 'hash');
end
