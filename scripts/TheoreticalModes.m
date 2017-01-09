% This script calculates the theoretical mode shapes at TAO mooring
% locations
function [] = TheoreticalModes()

  ticstart = tic;

  datadir = '../data/';
  woaTname = [datadir '/woa13_decav_t00_01v2.nc'];
  woaSname = [datadir '/woa13_decav_s00_01v2.nc'];
  woaMname = [datadir '/woa_landsea_01.msk'];
  etoponame = [datadir '/ETOPO2v2g_f4.nc4'];

  disp(' Loading ETOPO2v2 data.');
  etopo.x = ncread(etoponame, 'x');
  etopo.y = ncread(etoponame, 'y');
  etopo.z = ncread(etoponame, 'z');

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
  tao.lon = -1 * fliplr([95 110 125 140 155 170 180 ...
                      -165 -156 -147 -137]);
  tao.lat = [8 5 2 0 -2 -5 -8]; % +ve North

  nlon = length(tao.lon);
  nlat = length(tao.lat);

  %% Iterate!
  clear data
  for mm=1:nlon
      for nn=1:nlat

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

          % Locate (tao.lon,lat) on WOA grid and save
          ilat = find_approx(woa.Y, tao.lat(nn), 1);
          ilon = find_approx(woa.X, tao.lon(mm), 1);
          flatbot.lon(mm) = woa.X(ilon);
          flatbot.lat(nn) = woa.Y(ilat);

          fprintf(['TAO at (%.1f %s, %.1f %s) | ' ...
                   'WOA at (%.1f %s, %.1f %s)\n'], ...
                  abs(tao.lon(mm)), upper(lonstr), ...
                  abs(tao.lat(nn)), upper(latstr), ...
                  abs(flatbot.lon(mm)), upper(lonstr), ...
                  abs(flatbot.lat(nn)), upper(latstr));

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
              woa.indbot(ilon, ilat) = indbotT;
          else
              T = woa.temp(:,ilat,ilon); %squeeze(woa.temp(ilon, ilat, :));
              S = woa.sal(:,ilat,ilon); %squeeze(woa.sal(ilon, ilat, :));
          end

          % Get etopo depth at location
          etDepth = -AvgEtopoDepth(etopo, flatbot.lon(mm), ...
                                   flatbot.lat(nn));

          % find nearest depth level in WOA data
          indbot = find_approx(woa.Z, etDepth) + 1;
          if indbot > woa.indbot(ilon, ilat)
              warning(['WOA has bottom ' ...
                       num2str(woa.Z(indbot) - etDepth, '%.2f') ...
                       ' m shallower than etopo!']);
              indbot = woa.indbot(ilon, ilat);
          end

          Zmode = avg1(woa.Z);
          dtdz = -avg1(gradient(T, woa.Z)); %diff(T)./diff(woa.Z);
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

          flatbot.etopoDepth(mm,nn) = etDepth;
          flatbot.N2(mm,nn,:) = N2;
          flatbot.IdealWMode(mm,nn,:,:) = Vmode;
          flatbot.IdealTempMode(mm,nn,:,:) = Tmode;
          flatbot.Twoa{mm,nn} = T;
          flatbot.Swoa{mm,nn} = S;
          flatbot.Zwoa = woa.Z;
      end
  end

  flatbot.zTmode = Zmode;
  flatbot.comment = ['flatbot.IdealTempMode is the theoretical' ...
                     ' temperature mode on grid flatbot.zTmode. ' ...
                     'Similarly, flatbot.IdealWMode is the ' ...
                     'vertical velocity mode shape; ' ...
                     'flatbot.(T,S,Z)woa are WOA temp, salt, depth'];
  flatbot.hash = githash([mfilename('fullpath') '.m']);
  save('flat-bot-modes.mat', 'flatbot');
  toc(ticstart);
end

function [depth] = AvgEtopoDepth(etopo, lon, lat)

    ilon1 = find_approx(etopo.x, floor(lon));
    ilat1 = find_approx(etopo.y, floor(lat));

    ilon2 = find_approx(etopo.x, ceil(lon));
    ilat2 = find_approx(etopo.y, ceil(lat));

    Z = etopo.z(ilon1:ilon2, ilat1:ilat2);
    Z(Z >= -1) = NaN; % NaN out land

    depth = nanmean(Z(:));
end