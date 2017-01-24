% Get mode shape from data by band pass filtering around expected frequency
function [modes] = InferModeShape(opt, data, lonrange, latrange)

  ticstart = tic;

  disp('Inferring mode shapes');

  if ~exist('data', 'var'), error('Provide data.'); end
  if strcmpi(opt.name, 'montecarlo')
      isMonteCarlo = 1;
  else
      isMonteCarlo = 0;
  end

  nlon = length(data.lon);
  nlat = length(data.lat);

  if ~exist('lonrange', 'var')
      lonrange = 1:nlon;
  end
  if ~exist('latrange', 'var')
      latrange = 1:nlat;
  end

  modes.lon = data.lon;
  modes.lat = data.lat;
  modes.depth = data.depth;

  for mm=lonrange
      for nn=latrange
          if isempty(data.dht{mm,nn}) | isempty(data.T{mm,nn})
              continue;
          end

          dht = data.dht{mm,nn} - nanmean(data.dht{mm,nn});
          dht = fill_gap(dht, 'linear', opt.InterpGapLength);
          % First do Bandpass filtering
          dhtfilt = BandPass(dht, opt.filt);

          if opt.TagainstDHT
              [NoiseAmp, NoiseSlope] = EstimateNoiseSpectrum(dht, opt);
          end

          if opt.debug
              opt.hdbg = figure;
              hdbgax1 = subplot(211);
              PlotSpectrum(data.dht{mm,nn});
              PlotSpectrum(dhtfilt);
              linex(1./opt.filt.cutoff);
          end

          % iterate over depths and
          % 1. Filter temperature,
          % 2. Estimate background noise spectrum for temperature
          clear Tfilt
          if ~opt.TagainstDHT
              clear NoiseAmp NoiseSlope
          end

          for ii = 1:length(modes.depth{mm,nn})
              data.T{mm,nn}(ii,:) = fill_gap(data.T{mm,nn}(ii,:), ...
                                             'linear', opt.InterpGapLength);
              if opt.filter_temp
                  Tfilt(ii,:) = BandPass(data.T{mm,nn}(ii,:), opt.filt);
              else
                  Tfilt(ii,:) = data.T{mm,nn}(ii,:);
              end

              if ~opt.TagainstDHT
                  [NoiseAmp(ii), NoiseSlope(ii)] = ...
                      EstimateNoiseSpectrum(data.T{mm,nn}(ii,:), opt);
              end

              if opt.debug
                  axes(hdbgax1);
                  PlotSpectrum(data.T{mm,nn}(ii,:));
                  PlotSpectrum(Tfilt(ii,:));
              end
          end

          % Make sure I'm using same time interval for both variables
          % Needed because dynamic height is not available at all
          % time points with temperature measurement
          range = findCommonTimeRange(data.timedht{mm,nn}, ...
                                      data.timetemp{mm,nn});

          % Initialize to nans
          infer_mode = nan(size(modes.depth{mm,nn}));
          infer_mode_error = infer_mode;
          dof = infer_mode;
          corrcoeff = infer_mode;
          Tstd = infer_mode;

          Tstd = nanstd(Tfilt(:,range)')';
          [infer_mode, infer_mode_error, corrcoeff, ...
           intercept, stderr] ...
              = DoRegression(dhtfilt, Tfilt(:, range), ...
                             NoiseAmp, NoiseSlope, opt);

          % normalize mode shapes if not monte-carlo
          % first OLS
          if isMonteCarlo
              imnorm = 1;
          else
              imnorm = findNormalization(infer_mode(:,1));
          end
          modes.InferredModeOLS{mm,nn} = infer_mode(:,1)./imnorm;
          modes.InferredModeErrorOLS{mm,nn} = ...
              infer_mode_error(:,1)./imnorm;
          modes.StdErrorOLS{mm,nn} = stderr./imnorm;

          % now WTLS
          if isMonteCarlo
              imnorm = 1;
          else
              imnorm = findNormalization(infer_mode(:,2));
          end
          modes.InferredModeWTLS{mm,nn} = infer_mode(:,2)./imnorm;
          modes.InferredModeErrorWTLS{mm,nn} = ...
              infer_mode_error(:,2)./imnorm;

          if ~isMonteCarlo & any(cut_nan(modes.InferredModeOLS{mm,nn}) > 1)
              warning('OLS regression ampl > 1');
              keyboard;
          end

          % for low correlation, this might not work!
          % assert(all(cut_nan(modes.InferredModeWTLS{mm,nn}) <= 1), ...
          % 'WTLS regression ampl > 1');

          % figure out phase with depth
          [~,maxDepth] = max(modes.InferredModeOLS{mm,nn});
          phaselag = nan(length(modes.InferredModeOLS{mm,nn}), 1);
          omega = mean(opt.filt.cutoff/2);
          for ii=1:length(modes.InferredModeOLS{mm,nn})
              [c, lags, delay(ii)] = GappyCorrelation(Tfilt(ii,:), ...
                                                      Tfilt(maxDepth,:));

              % in degrees
              phaselag(ii) = atan2( tan(...
                  2*pi/omega * delay(ii)), 1) * 180/pi;
          end

          modes.phaselag{mm,nn} = phaselag;
          modes.intercept{mm,nn} = intercept;
          modes.corr{mm,nn} = corrcoeff;

          modes.timedht{mm,nn} = data.timedht;
          modes.dhtfilt{mm,nn} = dhtfilt;
          modes.Tstd{mm,nn} = Tstd;
      end
  end

  modes.comment = ['InferredMode(modes.lon,modes.lat,modes.depth)' ...
                   'is the ' ...
                   'mode structure inferred by regressing ' ...
                   'band passed dyn. ht against band passed' ...
                   'temperature. ' ...
                   'Tfilt = bandpassed temperature (cell array) | ' ...
                   'dhtfilt = band passed dynamic height | ' ...
                   'Tstd = standard dev of temp time series at depth'];

  % don't overwrite if *not* inferring at all locations
  if ~isMonteCarlo & ...
          isequal(lonrange, 1:nlon) & isequal(latrange, 1:nlat)
      hash = githash([mfilename('fullpath') '.m']);
      disp('Saving mode structures');
      save([opt.name '-' opt.filt.window '.mat'], ...
           'modes', 'opt', 'hash');
  end
  toc(ticstart);
end

function [range] = findCommonTimeRange(timedht, timetemp)

    start = find(timetemp == timedht(1));
    stop = find(timetemp == timedht(end));
    range = start:stop;

    assert(all(timetemp(range) == timedht));
end

function [norm] = findNormalization(shape)

    [norm,imax] = nanmax(shape);

    % if imax > 17
    %     % max is too deep, generally happens with WTLS
    %     [norm,imax] = nanmax(shape(1:imax-1));
    % end

    % Normalize so that mode shape is 1 at water depth
    % closest to where the theoretical mode maximum is.
    % [~, izmaxwoa] = max(flatbot.IdealTempMode(mm,nn,:, ...
    %                     opt.nmode));
    % zind = find_approx(modes.depth{mm,nn}, ...
    %                    flatbot.Zwoa(izmaxwoa), 5);

    % % normalize but account for NaNs and 0s
    % imnorm = cut_nan(fillnan(infer_mode(zind, 1), 0));
    % imnorm = max(imnorm);
end