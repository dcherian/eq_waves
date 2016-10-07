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

  ticstart = tic;
  datadir = '../data/';
  woaTname = [datadir '/woa13_decav_t00_01v2.nc'];
  woaSname = [datadir '/woa13_decav_s00_01v2.nc'];
  woaMname = [datadir '/woa_landsea_01.msk'];
  etoponame = [datadir '/ETOPO2v2g_f4.nc4'];

  fprintf('\n Loading WOA 13 data \n\n');
  woa.X = double(ncread(woaSname, 'lon'));
  woa.Y = double(ncread(woaSname, 'lat'));
  woa.Z = double(ncread(woaSname, 'depth'));
  woa.temp = permute(double(ncread(woaTname, 't_an')), [3 2 1]);
  woa.sal = permute(double(ncread(woaSname, 's_an')), [3 2 1]);

  % woaMname contains mask with index of standard depth level that
  % is closest to the bottom. Map this.
  woamask = csvread(woaMname, 2);
  woa.indbot = griddata(woamask(:,2), woamask(:,1), woamask(:,3), ...
                        woa.X', woa.Y, 'nearest')';

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

          % TAO filenames
          fnamet = [datadir 'temp/t',   num2str(abs(modes.lat(nn))), ...
                    latstr,num2str(abs(modes.lon(mm))),lonstr,'_dy.cdf'];
          fnameh = [datadir 'dynht/dyn',num2str(abs(modes.lat(nn))), ...
                    latstr,num2str(abs(modes.lon(mm))),lonstr,'_dy.cdf'];

          % make sure observations exist before continuing
          if ~exist(fnamet,'file'), continue; end
          if ~exist(fnameh,'file'), continue; end

          %% Step 1. Calculate theoretical modes

          % Locate (modes.lon,lat)
          ilat = find_approx(woa.Y, modes.lat(nn), 1);
          ilon = find_approx(woa.X, modes.lon(mm), 1);

          fprintf(['TAO at (%.1f %s, %.1f %s) | ' ...
                   'WOA at (%.1f %s, %.1f %s)\n'], ...
                  abs(modes.lon(mm)), upper(lonstr), ...
                  abs(modes.lat(nn)), upper(latstr), ...
                  abs(woa.X(ilon)), upper(lonstr), ...
                  abs(woa.Y(ilat)), upper(latstr));

          if all(isnan(woa.temp(:,ilat,ilon)))
              warning('No climatological data! Averaging...');
              T = nanmean([woa.temp(:,ilat-1,ilon) ...
                           woa.temp(:,ilat+1,ilon) ...
                           woa.temp(:,ilat,ilon-1) ...
                           woa.temp(:,ilat,ilon+1)], 2);
              S = nanmean([woa.sal(:,ilat-1,ilon) ...
                           woa.sal(:,ilat+1,ilon) ...
                           woa.sal(:,ilat,ilon-1) ...
                           woa.sal(:,ilat,ilon+1)], 2);
              indbotT = find(isnan(T) == 1, 1, 'first');
              indbotS = find(isnan(S) == 1, 1, 'first');
              assert(indbotT == indbotS);
              indbot = indbotT;
          else
              T = woa.temp(:,ilat,ilon); %squeeze(woa.temp(ilon, ilat, :));
              S = woa.sal(:,ilat,ilon); %squeeze(woa.sal(ilon, ilat, :));
              % figure out water depth from WOA land-sea mask
              indbot = woa.indbot(ilon, ilat);
          end

          if indbot > length(T) % all valid data
              indbot = length(T)+1;
          else
              % make sure I have correct bottom
              % and that griddata has not screwed up
              assert(isnan(T(indbot)) & isnan(S(indbot)));
          end

          Zmode = avg1(woa.Z);
          dtdz = avg1(gradient(T, woa.Z)); %diff(T)./diff(woa.Z);
          dtdz(indbot-1:end) = NaN;

          N2 = bfrq(S,T,woa.Z,modes.lat(nn));%10^(-6)*ones(32,1);
          N2(N2 < 0) = min(abs(N2(:)));
          N2(indbot-1:end) = NaN;

          % Vmode = vertical velocity mode shape
          [Vmode, ~, ~] = vertmode(N2(1:indbot-2), ...
                                   woa.Z(1:indbot-1),3,0);
          Vmode(indbot-1:length(Zmode),:) = NaN; % extend to bottom
          % calculate temperature mode shape
          Tmode = Vmode .* repmat(dtdz,1,size(Vmode,2));

          % make sure maximum is always positive.
          sgn = -1 * ~(max(Tmode) == max(abs(Tmode)));
          sgn(sgn == 0) = 1;
          Tmode = Tmode .* sgn;

          %% Filter out TAO data

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
              data.Tfilt{mm,nn,ii,:} = Tfilt(ii,:);
              Tstd(ii) = nanstd(Tfilt(ii,:));

              % find all nan's in both datasets
              Treduced = Tfilt(ii, range);
              mask = ~(isnan(dhtfilt) | isnan(Treduced));

              Treduced(~mask) = NaN;
              dhtfiltreg = dhtfilt; % copy for regression
              dhtfiltreg(~mask) = NaN;

              if all(mask == 0)
                  infer_mode(ii) = NaN;
                  infer_mode_error(ii,1) = NaN;
                  continue;
              end

              % regress to find mode shape
              % infer_mode(ii) = dhtfilt(mask)' \ Treduced(mask)';
              % [infer_mode(ii), bint] = ...
              %     regress(Treduced(mask)', dhtfilt(mask)');
              % infer_mode_error(ii,1) = bint(2) - infer_mode(ii);

              % [coef, conf, dof(ii)] = dcregress(dhtfiltreg', Treduced', ...
              %                                   [], 0);

              rr = mf_wtls(dhtfiltreg', Treduced', nanstd(dhtfiltreg), ...
                           nanstd(Treduced), 0);
              coef(1) = rr(3);
              coef(2) = rr(1);
              conf(1) = rr(4);
              conf(2) = rr(2);

              infer_mode(ii) = coef(2);
              infer_mode_error(ii) = conf(2);
              corrcoeff(ii) = min(min( ...
                  corrcoef(dhtfilt(mask)', Treduced(mask)')));

              if abs(corrcoeff(ii)) <= corr_sig(dof(ii), 0.95)
                  % 0 means insignificant, NaN means no data.
                  corrcoeff(ii) = 0;
              end

              if isnan(dof(ii))
                  dof(ii) = min([calcdof(dhtfiltreg) ...
                                 calcdof(Treduced)]);
              end

              if length(cut_nan(dhtfiltreg(1:Nsamp:end))) <= 2 | ...
                      length(cut_nan(isnan(Treduced(1:Nsamp:end)))) <= 2
                  continue;
              end

              rr = mf_wtls(dhtfiltreg(1:Nsamp:end)', ...
                           Treduced(1:Nsamp:end)', ...
                           nanstd(dhtfiltreg), nanstd(Treduced), 0);
              coef(1) = rr(3);
              coef(2) = rr(1);
              conf(1) = rr(4);
              conf(2) = rr(2);

              infer_mode(ii) = coef(2);
              infer_mode_error(ii) = conf(2);

>>>>>>> 8112abf... reorder corr sig.
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

          modes.InferredMode{mm,nn} = infer_mode./nanmax(infer_mode);
          modes.InferredModeError{mm,nn} = infer_mode_error./nanmax(abs(infer_mode));
          modes.IdealTempMode(mm,nn,:,:) = Tmode;
          modes.dof{mm,nn} = dof;
          modes.corr{mm,nn} = corrcoeff;

          data.timedht{mm,nn} = tao.timedht;
          data.dhtfilt{mm,nn} = dhtfilt;
          data.Tstd{mm,nn} = Tstd;
          data.Twoa{mm,nn} = T;
          data.Swoa{mm,nn} = S;
          data.Zwoa = woa.Z;
      end
  end

  modes.zTmode = Zmode;

  hash = githash([mfilename('fullpath') '.m']);
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

  save([opt.name '-' opt.filt.window '.mat'], 'modes', 'data', 'opt', 'hash');
  toc(ticstart);
end
