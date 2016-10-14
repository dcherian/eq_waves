% This script calculates the theoretical mode shapes at TAO mooring
% locations
function [] = TheoreticalModes()

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

  % edited lon values to make neater subplots ...
  % (some lons dont have enough data)
  % values in +ve East.
  flatbot.lon = -1 * fliplr([95 110 125 140 155 170 180 ...
                      -165 -156 -147 -137]);
  flatbot.lat = [8 5 2 0 -2 -5 -8]; % +ve North

  nlon = length(flatbot.lon);
  nlat = length(flatbot.lat);

  %% Iterate!
  clear data
  for mm=1:nlon
      for nn=1:nlat

          if flatbot.lon(mm) > 0
              lonstr = 'e';
          else
              lonstr = 'w';
          end

          if flatbot.lat(nn) < 0
              latstr = 's';
          else
              latstr = 'n';
          end

          % TAO filenames
          fnamet = [datadir 'temp/t',   num2str(abs(flatbot.lat(nn))), ...
                    latstr,num2str(abs(flatbot.lon(mm))),lonstr,'_dy.cdf'];
          fnameh = [datadir 'dynht/dyn',num2str(abs(flatbot.lat(nn))), ...
                    latstr,num2str(abs(flatbot.lon(mm))),lonstr,'_dy.cdf'];

          % make sure observations exist before continuing
          if ~exist(fnamet,'file'), continue; end
          if ~exist(fnameh,'file'), continue; end

          % Locate (flatbot.lon,lat)
          ilat = find_approx(woa.Y, flatbot.lat(nn), 1);
          ilon = find_approx(woa.X, flatbot.lon(mm), 1);

          fprintf(['TAO at (%.1f %s, %.1f %s) | ' ...
                   'WOA at (%.1f %s, %.1f %s)\n'], ...
                  abs(flatbot.lon(mm)), upper(lonstr), ...
                  abs(flatbot.lat(nn)), upper(latstr), ...
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

          N2 = bfrq(S,T,woa.Z,flatbot.lat(nn));%10^(-6)*ones(32,1);
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
          Tmode = bsxfun(@times, Tmode, sgn);

          flatbot.IdealTempMode(mm,nn,:,:) = Tmode;
          flatbot.Twoa{mm,nn} = T;
          flatbot.Swoa{mm,nn} = S;
          flatbot.Zwoa = woa.Z;
      end
  end

  flatbot.zTmode = Zmode;
  flatbot.comment = ['flatbot.IdealTempMode is the theoretical' ...
                     ' temperature mode on grid flatbot.zTmode. ' ...
                     'flatbot(T,S,Z)woa are WOA temp, salt, depth'];
  flatbot.hash = githash([mfilename('fullpath') '.m']);
  save('flat-bot-modes.mat', 'flatbot');
  toc(ticstart);
end