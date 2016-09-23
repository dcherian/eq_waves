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


  fprintf('\n Loading WOA 13 data \n\n');
  woa.X = double(ncread([datadir '/woa13_decav_s00_01v2.nc'], 'lon'));
  woa.Y = double(ncread([datadir '/woa13_decav_s00_01v2.nc'], 'lat'));
  woa.Z = double(ncread([datadir '/woa13_decav_s00_01v2.nc'], 'depth'));
  woa.temp = permute(double(ncread([datadir '/woa13_decav_t00_01v2.nc'], ...
                                   't_an')), [3 2 1]);
  woa.sal = permute(double(ncread([datadir '/woa13_decav_s00_01v2.nc'], ...
                                  's_an')), [3 2 1]);

  modes.lat = [8 5 2 0 -2 -5 -8]; % N
  % modes.lon = [95 110 125 140 137 147 155 156 165 170 180]; % W
  % edited lon values to make neater subplots ...
  % (some lons dont have enough data)
  % values in +ve East.
  modes.lon = -1 * fliplr([95 110 125 140 155 170 180 -165 -156 -147 -137]);

  nlon = length(modes.lon);
  nlat = length(modes.lat);

  %% Iterate!
  clear data
  for mm=1:nlon
      for nn=1:nlat

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

          %% Step 1. Calculate theoretical modes

          % Locate (modes.lon,lat)
          ilat = find_approx(woa.Y, modes.lat(nn), 1);
          ilon = find_approx(woa.X, modes.lon(mm), 1);

          fprintf('\n Actual location = (%.1f %s, %.1f %s) \n', ...
                  abs(woa.X(ilon)), upper(lonstr), ...
                  abs(woa.Y(ilat)), upper(latstr));

          Zmode = avg1(woa.Z);
          T = woa.temp(:,ilat,ilon); %squeeze(woa.temp(ilon, ilat, :));
          S = woa.sal(:,ilat,ilon); %squeeze(woa.sal(ilon, ilat, :));
          if all(isnan(T)) | all(isnan(S))
              warning('no climatological data!');
              continue;
          end
          dtdz = diff(T)./diff(woa.Z);

          N2 = bfrq(S,T,woa.Z,modes.lat(nn));%10^(-6)*ones(32,1);
          N2(N2 < 0) = min(abs(N2(:)));

          [Vmode, Hmode, c] = vertmode(N2,woa.Z,3,0);
          % calculate temperature mode shape
          Tmode = Vmode .* repmat(dtdz,1,size(Vmode,2));

          % make sure maximum is always positive.
          sgn = -1 * ~(max(Tmode) == max(abs(Tmode)));
          sgn(sgn == 0) = 1;
          Tmode = Tmode .* sgn;

          %% Filter out TAO data

          % clear tbuoy depth dht Tfilt Tfilt1 dhtfilt

          % read in TAO measurements
          fnamet = [datadir 'temp/t',   num2str(abs(modes.lat(nn))), ...
                    latstr,num2str(abs(modes.lon(mm))),lonstr,'_dy.cdf'];
          fnameh = [datadir 'dynht/dyn',num2str(abs(modes.lat(nn))), ...
                    latstr,num2str(abs(modes.lon(mm))),lonstr,'_dy.cdf'];

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
          dhtfilt = BandPass(tao.dht, opt.filt);

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
          infer_mode_error = nan(size(modes.depth{mm,nn}));
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
              data.Tfilt{mm,nn,ii,:} = Tfilt(ii,:);
              Tstd(ii) = nanstd(Tfilt(ii,:));

              % find all nan's in both datasets
              Treduced = Tfilt(ii, range);
              mask = ~(isnan(dhtfilt) | isnan(Treduced));

              if all(mask == 0)
                  infer_mode(ii) = NaN;
                  infer_mode_error(ii,1) = NaN;
                  continue;
              end

              % regress to find mode shape
              infer_mode(ii) = dhtfilt(mask)' \ Treduced(mask)';
              [infer_mode(ii), bint] = ...
                  regress(Treduced(mask)', dhtfilt(mask)');
              infer_mode_error(ii,1) = bint(2) - infer_mode(ii);

              if opt.debug
                  figure(hdbg);
                  hdbgax2 = subplot(212);
                  plot(dhtfilt(mask)); hold on;
                  plot(Treduced(mask));
                  title('filtered time series for regression');
                  keyboard;
              end
              %Imode = fill_gap(dhtfilt(range)','linear',15)\fill_gap(Tfilt(:,range)','linear',15);
          end
          modes.InferredMode{mm,nn} = infer_mode./nanmax(abs(infer_mode));
          modes.InferredModeError{mm,nn} = infer_mode_error./nanmax(abs(infer_mode));
          modes.IdealTempMode(mm,nn,:,:) = Tmode;
          data.Tstd{mm,nn} = Tstd;

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
