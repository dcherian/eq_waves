% Get mode shape from data by band pass filtering around expected frequency
% Deepak Cherian

% All below fixed - 25th July 2012
% Issues (30th Sept. 2011)
% 1) Present filtering method results in crap in the first and last two points.
% 2) The inferred mode isnt being saved anywhere.
% 3) The main issue is finding a suitable time interval to run the
% regression in so that we get max. use of available data. Need to remember
% how the old TAO data scripts work.

%% Data & Parameters

clr
load 'data\dynht.mat'

if exist('temp','var') ~= 1 && exist('sal','var') ~=1
    fprintf('\n Loading WOA 05 data \n\n');
    load 'data\woa05.mat';
end

modes.lat = [8 5 2 0 -2 -5 -8]; % N
%modes.lon = [95 110 125 140 137 147 155 156 165 170 180]; % W
% edited lon values to make neater subplots (some lons dont have enough data)
modes.lon = fliplr([95 110 125 140 155 170 180]);

modes.depth = [1 25 50 75 100 125 150 200 250 300 500]'; % standard depths
n_modes = 3;
fontSize = [12 12 14];

nlon = length(modes.lon);
nlat = length(modes.lat);

%% Iterate!

for mm=1:nlon
    for nn=1:nlat
        %% Calculate theoretical modes

        % Locate (modes.lon,lat)     
        ilon = find_approx(X,360-modes.lon(mm),1);
        ilat = find_approx(Y,modes.lat(nn),1);

        % Verify location
        fprintf('\n Actual location = (%.1f W,%.1f N) \n', 360-X(ilon), Y(ilat));

        % Interpolate to denser grid.
        Zmode = avg1(Z);
        T = temp(:,ilat,ilon);
        S = sal(:,ilat,ilon);
        dtdz = diff(T)./diff(Z);
        %P = sw_pres(modes.depth,modes.lat);

        N2 = bfrq(S,T,Z,modes.lat(nn));%10^(-6)*ones(32,1);
        [Vmode, Hmode, c] = vertmode(N2,Z,n_modes,0);
        Tmode = Vmode .* repmat(dtdz,1,size(Vmode,2));

        %% Filter

        clear vars atts dims tbuoy depth dht tavg tavg1 dhtavg

        if modes.lat(nn) < 0
            fnamet = ['data\temp\t',   num2str(abs(modes.lat(nn))),'s',num2str(modes.lon(mm)),'w_dy.cdf'];
            fnameh = ['data\dynht\dyn',num2str(abs(modes.lat(nn))),'s',num2str(modes.lon(mm)),'w_dy.cdf'];
        else
            fnamet = ['data\temp\t',   num2str(abs(modes.lat(nn))),'n',num2str(modes.lon(mm)),'w_dy.cdf'];            
            fnameh = ['data\dynht\dyn',num2str(abs(modes.lat(nn))),'n',num2str(modes.lon(mm)),'w_dy.cdf'];
        end
        
        % Read & interpolate temp
        if ~exist(fnamet,'file'), continue; end
        %[vars atts dims] = ncdfread(fnamet);
        tbuoy = addnan(squeeze(ncread(fnamet,'T_20')),100);
        depth = squeeze(ncread(fnamet,'depth'));
        tbuoy = interp1(depth,tbuoy,modes.depth);
        
        % Read dynamic ht
        if ~exist(fnameh,'file'), continue; end
        %clear vars atts dims
        %[vars atts dims] = ncdfread(fnameh);
        dht = addnan(squeeze(ncread(fnameh,'DYN_13')),1000)';

        % Running averages 6 day , 12 day and subtract
        windows = [6 12];
        dhtavg = conv_band_pass(dht,windows);
        
        % Plot 'mode' structure
        clear Imode
        n_mode = 2;
        range = 1:length(dhtavg);
        
        % iterate over standard depths
        for ii = 1:11
            % band pass temperature data
            tavg(ii,:) = conv_band_pass(tbuoy(ii,:),windows);
            
            % find all nan's in both datasets, negate that mask -> NaN locations are 0
            % replace 0 with NaN and multiple data and remove those
            mask = fillnan(double(~(isnan(dhtavg(range)) | isnan(tavg(ii,range)))),0);
            modes.Imode(mm,nn,ii) = cut_nan(dhtavg(range).*mask)'\cut_nan(tavg(ii,range).*mask)';
        end
        %Imode = fill_gap(dhtavg(range)','linear',15)\fill_gap(tavg(:,range)','linear',15);
        ind500 = find_approx(Zmode,500,1);
        modes.Imode(mm,nn,:) = modes.Imode(mm,nn,:)./max(modes.Imode(mm,nn,:));
        modes.Tmode(mm,nn,:) = Tmode(:,n_mode);
        % now plot
%         figure(1); clf
%         plot(squeeze(modes.Imode(mm,nn,:)),modes.depth,'bo-');
%         hold on
%         plot(abs(Tmode(:,n_mode)')./max(abs(Tmode(1:ind500,n_mode))),Zmode,'r'); % 
%         title(['(',num2str(modes.lon(mm)), 'W, ', num2str(modes.lat(nn)) , 'N)']);
%         %legendflex({'Inferred mode', 'Theoretical mode'},'anchor',{'s','w'},'nrow',1,'buffer',[40 40],'fontsize',12);
%         legend('Inferred mode','Theoretical mode','Location','SouthEast');
%         plot(zeros(size(modes.depth)),modes.depth,'k-')
%         %xlim([-0.2 1.2]);
%         ylim([modes.depth(1) modes.depth(end)]);
%         revz;
%         beautify(fontSize);
%         
%         printname = ['images\', num2str(modes.lon(mm)), 'W', num2str(modes.lat(nn)) , 'N', '.png'];
%         export_fig('-nocrop',printname);
    end
end

%% combined plot
figure(2); 
set(gcf,'Position',[0 0 1600 900]);
xlimits = [-0.2 1.2];
fontSize = [11 13 14];

for mm=1:nlon
    for nn=1:nlat 
        if (modes.Imode(mm,nn,:) == 0), continue; end
        subplot_index = sub2ind([nlon nlat],mm,nn);
        subplot(nlat,nlon,subplot_index);
        plot(squeeze(modes.Imode(mm,nn,:)),modes.depth,'bo-');
        hold on
        plot(abs(Tmode(:,n_mode)')./max(abs(Tmode(1:ind500,n_mode))),Zmode,'r'); % 
        plot(zeros(size(modes.depth)),modes.depth,'-','Color',[0.3 0.3 0.3])
        xlim(xlimits);
        ylim([modes.depth(1) modes.depth(end)]);
        set(gca,'YTick',[modes.depth(1) get(gca,'YTick')]);
        revz;
        
        if mod(subplot_index,nlon) == 1
            ylabel([num2str(modes.lat(nn))  'N']); 
        else
            set(gca,'YTickLabel',[]);
        end
        if (subplot_index >= nlon*nlat - nlon+1)
            xlabel([num2str(modes.lon(mm)) 'W']); 
        else
            set(gca,'XTickLabel',[]);
        end
        beautify(fontSize);
        set(gca,'TickLength',[0.03 0.03]);
        box off
    end
end

% fix missing lon labels
missing_lon = 140;
figure(2)
for ii = 1:length(missing_lon)
    subplot(nlat,nlon,sub2ind([nlon nlat],find(modes.lon == missing_lon(ii)),nlat))
    xlabel('140 W');
    set(gca,'YTickLabel',[]);
    beautify(fontSize);
    set(gca,'YColor',[1 1 1]);
    box off
end

[ax,h] = suplabel('Red = Theoretical Mode, Blue = Inferred Mode','t');
set(h,'FontSize',fontSize(3),'FontName','AvantGarde');

export_fig('-nocrop','-r150','combined.png');

%% save to file
modes.comment = ['Imode(modes.lon,modes.lat,modes.depth) is the mode structure inferred by regressing band passed (6 day - 12 day) dyn. ht against band passed temperature.' ...
                    ' modes.Tmode is the theoretical temperature mode']; 
save Imode.mat modes
