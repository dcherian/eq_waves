% Get mode shape from data by band pass filtering around expected frequency
function [modes] = InferModeShape(opt, data, lonrange, latrange)

  ticstart = tic;

  disp('Inferring mode shapes');

  if ~exist('data', 'var'), error('Provide data.'); end

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

  if strcmpi(data.name, 'tao')
      if exist('bounds.mat', 'file')
          load bounds.mat
      else
          mbound = [];
      end
  end

  for mm=lonrange
      for nn=latrange
          disp([mm nn])

          if isempty(data.dht{mm,nn}) | isempty(data.T{mm,nn})
              continue;
          end

          dht = data.dht{mm,nn} - nanmean(data.dht{mm,nn});
          dht = fill_gap(dht, 'linear', opt.InterpGapLength);
          % First do Bandpass filtering
          dhtfilt = BandPass(dht, opt.filt);

          if ~opt.TagainstDHT
              % use longest valid subset
              [noise.amp, noise.slope] = ...
                  EstimateNoiseSpectrum(dht, opt, 0, [], []);
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

          for ii = 1:length(modes.depth{mm,nn})
              data.T{mm,nn}(ii,:) = fill_gap(data.T{mm,nn}(ii,:), ...
                                             'linear', opt.InterpGapLength);
              if opt.filter_temp
                  Tfilt(ii,:) = BandPass(data.T{mm,nn}(ii,:), opt.filt);
              else
                  Tfilt(ii,:) = data.T{mm,nn}(ii,:);
              end

              if opt.TagainstDHT
                  [noise.amp(ii), noise.slope(ii)] = ...
                      EstimateNoiseSpectrum(data.T{mm,nn}(ii,:), ...
                                            opt, 0, [], []);
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

          sigslope.m = mbound{mm,nn};
          sigslope.mdist = mdist{mm,nn};

          Tstd = nanstd(Tfilt(:,range)')';
          [infer_mode, infer_mode_error, corrcoeff, ...
           intercept, stderr, dof] ...
              = DoRegression(dhtfilt, Tfilt(:, range), ...
                             noise, sigslope, opt, dht);

          % normalize mode shapes if not monte-carlo
          % first OLS
          imnorm = findNormalization(infer_mode(:,1));
          modes.InferredModeOLS{mm,nn} = infer_mode(:,1)./imnorm;
          modes.OLSnorm(mm,nn) = imnorm;
          modes.InferredModeErrorOLS{mm,nn} = ...
              infer_mode_error(:,1)./imnorm;
          modes.StdErrorOLS{mm,nn} = stderr./imnorm;
          modes.dof{mm,nn} = dof;

          % now WTLS
          imnorm = findNormalization(infer_mode(:,2));
          modes.InferredModeWTLS{mm,nn} = infer_mode(:,2)./imnorm;
          modes.InferredModeErrorWTLS{mm,nn} = ...
              infer_mode_error(:,2)./imnorm;

          if any(cut_nan(modes.InferredModeOLS{mm,nn}) > 1)
              warning('OLS regression ampl > 1');
              keyboard;
          end

          % for low correlation, this might not work!
          % assert(all(cut_nan(modes.InferredModeWTLS{mm,nn}) <= 1), ...
          % 'WTLS regression ampl > 1');

          % figure out phase with depth
          [~,maxDepth] = max(modes.InferredModeOLS{mm,nn});
          phaselag = nan(length(modes.InferredModeOLS{mm,nn}), 1);
          omega = 2*pi/mean(opt.filt.cutoff/2);
          for ii=1:length(modes.InferredModeOLS{mm,nn})
              [c, lags, delay(ii)] = GappyCorrelation(Tfilt(ii,:), ...
                                                      Tfilt(maxDepth,:));

              % phaselag in degrees; -180 <= phaselag <= 180
              % delay is in days; omega is in cpd.
              phaselag(ii) = atan2(tan(omega * delay(ii)), 1) * 180/pi;
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
  if isequal(lonrange, 1:nlon) & isequal(latrange, 1:nlat)
      hash = githash([mfilename('fullpath') '.m']);
      disp('Saving mode structures');
      save([opt.name '-' opt.filt.window '.mat'], ...
           'modes', 'opt', 'hash');
  end
  toc(ticstart);
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