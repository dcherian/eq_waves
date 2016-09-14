% Get mode shape from data by band pass filtering around expected frequency
% Deepak Cherian

% All below fixed - 25th July 2012
% Issues (30th Sept. 2011)
% 1) Present filtering method results in crap in the first and last two points.
% 2) The inferred mode isnt being saved anywhere.
% 3) The main issue is finding a suitable time interval to run the
% regression in so that we get max. use of available data. Need to remember
% how the old TAO data scripts work.

% This is how the script works for each point in the TAO array
% 1) Load WOA data and calculate theoretical 0 mean vel. + flat bottom mode
% shapes. The T mode is obtained by multiplying the vertical mode by dT/dz
% (from WOA).
% 2) Read TAO temp data, ***interpolate*** to standard depth grid.
% 3) Read dynamic height data and band pass filter around given band (now 6 
% and 14 days)
% 4) At each depth, band pass filter interpolated temp data, find missing 
% data in both DHT & temp, and remove them. Then, do the inversion and
% calculate "inferred mode" (Imode) again at each depth.
% 5) Normalize by *max amplitude* and save data.

%% Data & Parameters

clear

% main options
butterworth = 1; % butterworth filter if 1, else running mean
filter_temp = 1; % filter temperature also?
windows = [6 12]; % (days) band pass filter windows
n_modes = 3; % number of modes to calculate
n_mode = 2; % which theoretical mode am I looking for?

debug = 0; % debugging spectrum plots

datadir = '../data/';
load([datadir 'dynht.mat']);

if exist('temp','var') ~= 1 && exist('sal','var') ~=1
    fprintf('\n Loading WOA 05 data \n\n');
    load '../data/woa05.mat';
end

modes.lat = [8 5 2 0 -2 -5 -8]; % N
%modes.lon = [95 110 125 140 137 147 155 156 165 170 180]; % W
% edited lon values to make neater subplots (some lons dont have enough data)
modes.lon = fliplr([95 110 125 140 155 170 180]);

%modes.depth = [1 25 50 75 100 125 150 200 250 300 500]'; % standard depths
fontSize = [12 12 14];

nlon = length(modes.lon);
nlat = length(modes.lat);

%% Iterate!
clear data
for mm=1:nlon
    for nn=1:nlat
        %% Calculate theoretical modes

        % Locate (modes.lon,lat)     
        ilon = find_approx(X,360-modes.lon(mm),1);
        ilat = find_approx(Y,modes.lat(nn),1);

        % Verify location
        fprintf('\n Actual location = (%.1f W,%.1f N) \n', 360-X(ilon), Y(ilat));

        % Interpolate to denser grid.
        Zmode = avg1(Z); % Z from woa05.mat
        T = temp(:,ilat,ilon);
        S = sal(:,ilat,ilon);
        dtdz = diff(T)./diff(Z);
        %P = sw_pres(modes.depth,modes.lat);

        N2 = bfrq(S,T,Z,modes.lat(nn));%10^(-6)*ones(32,1);
        [Vmode, Hmode, c] = vertmode(N2,Z,n_modes,0);
        % calculate temperature mode shape
        Tmode = Vmode .* repmat(dtdz,1,size(Vmode,2));

        %% Filter

        clear vars atts dims tbuoy depth dht tavg tavg1 dhtavg

        if modes.lat(nn) < 0
            fnamet = [datadir 'temp/t',   num2str(abs(modes.lat(nn))),'s',num2str(modes.lon(mm)),'w_dy.cdf'];
            fnameh = [datadir 'dynht/dyn',num2str(abs(modes.lat(nn))),'s',num2str(modes.lon(mm)),'w_dy.cdf'];
        else
            fnamet = [datadir 'temp/t',   num2str(abs(modes.lat(nn))),'n',num2str(modes.lon(mm)),'w_dy.cdf'];            
            fnameh = [datadir 'dynht/dyn',num2str(abs(modes.lat(nn))),'n',num2str(modes.lon(mm)),'w_dy.cdf'];
        end
        
        % Read & interpolate temp
        if ~exist(fnamet,'file'), continue; end
        %[vars atts dims] = ncdfread(fnamet);
        tbuoy = addnan(squeeze(ncread(fnamet,'T_20')),100);
        modes.depth{mm,nn} = squeeze(ncread(fnamet,'depth'));
        data.depth{mm,nn}  = modes.depth{mm,nn};
        
        % interpolate buoy profiles onto standard depth - COULD BE IMPROVED
        %tbuoy = interp1(depth,tbuoy,modes.depth);
        
        % Read dynamic ht
        if ~exist(fnameh,'file'), continue; end
        %clear vars atts dims
        %[vars atts dims] = ncdfread(fnameh);

        dht = double(addnan(squeeze(ncread(fnameh,'DYN_13')),1000))';

        dht = dht - nanmean(dht);

        if butterworth
            [b,a] = butter(12, sort(2*pi./windows/(2*pi/2)), 'bandpass');
            dhtavg = filter(b,a,dht);
        else
            % Running averages windows(1) , windows(2) day and subtract
            dhtavg = conv_band_pass(dht,windows);
        end

        if debug
            PlotSpectrum(dht);
            PlotSpectrum(dhtavg);
            linex(1./windows);
            keyboard;
        end

        % Plot 'mode' structure        
        %        range = 1:length(dhtavg);
        
        infer_mode = nan(size(modes.depth{mm,nn}));
        tstd = infer_mode;
        % iterate over standard depths
        for ii = 1:length(modes.depth{mm,nn})
            if filter_temp
                % band pass temperature data
                if butterworth
                    tavg(ii,:) = filter(b,a,tbuoy(ii,:));
                else
                    tavg(ii,:) = conv_band_pass(tbuoy(ii,:),windows);
                end
            else
                tavg(ii,:) = tbuoy(ii,:);
            end

            if debug
                PlotSpectrum(tavg(ii,:));
                keyboard;
            end

            % save temperature series for reference and calculate std dev 
            data.tavg{mm,nn,ii,:} = tavg(ii,:);
            tstd(ii) = nanstd(tavg(ii,:));
            
            % find all nan's in both datasets, negate that mask -> NaN locations are 0
            % replace 0 with NaN and multiple data and remove those
            mask = fillnan(double(~(isnan(dhtavg(range)) | isnan(tavg(ii,range)))),0);
            infer_mode(ii) = cut_nan(dhtavg(range).*mask)' ...
                \ cut_nan(tavg(ii,range).*mask)';
        end
        %Imode = fill_gap(dhtavg(range)','linear',15)\fill_gap(tavg(:,range)','linear',15);
        
        modes.Imode{mm,nn} = infer_mode./max(infer_mode);
        modes.Tmode(mm,nn,:) = Tmode(:,n_mode);
        data.tstd{mm,nn} = tstd;
        
        % now plot - moved to plot_modes.m - this script now just creates a
        % mat file that contains data to be plotted
    end
end

modes.zTmode = Zmode;

data.hash = githash([mfilename('fullpath') '.m']);
data.dhtavg = dhtavg;
data.comment = ['tavg = bandpassed temperature (cell array) | ' ...
                'dhtavg = band passed dynamic height' ...
                'tstd = standard dev of temp time series at depth'];

%% save to file
modes.comment = ['Imode(modes.lon,modes.lat,modes.depth) is the mode structure inferred by regressing band passed (6 day - 12 day) dyn. ht against band passed temperature.' ...
                    ' modes.Tmode is the theoretical temperature mode on grid modes.zTmode']; 
save Imode.mat modes data
