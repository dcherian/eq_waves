% Uses Combined Data
% Date: 4th jan 2010

% compare mean T,S

%% Prepare data

n_modes = 30; % number of BAROCLINIC modes

lat = 0; % N
lon = 360-23; % E, APPROXIMATE

if exist('U_all','var') ~= 1 
     fprintf('\n Loading combined data ...\n');
     load 'E:\Work\Summer 2010\eql_Atl_combined_687_842_redone_35.mat';
end


%% Get normalized vertical modes for location (lon,lat)

if exist('temp','var') ~= 1 && exist('sal','var') ~=1
    fprintf('\n Loading WOA 05 data');
    load woa05.mat;
end

% Locate (lon,lat) % ASSUMES EITHER INTEGER OR ___.5 N/E
if mod(lon,1) == 0
    x = find(X < lon+1 & X > lon-1);
else
    x = find(X == lon);
end
y = find(Y < lat+1 & Y > lat-1);
ilon = x(1); ilat = y(1); % locations in WOA grid
fprintf('\n Actual location = (%.1f W,%.1f N)', 360-X(ilon), Y(ilat));

% Interpolate to denser grid.
T1 = temp(:,ilat,ilon);
S1 = sal(:,ilat,ilon);

% Use first missing value depth in climatological data to fix bottom
% Needed because vertmode.m calculates modes for full water column
miss = find(isnan(T1) == 1);
if isempty(miss), miss = length(Z); end

% Zext -> input grid extended to surface and bottom (bottom as in WOA'05)
% depths -> grid such that midpoints are at Zext
% Zin = z1;%sw_dpth(pgrid,lat)';
dP = 4;% pgrid(2)-pgrid(1); % grid spacing in dbar
Zext = sw_dpth(dP:dP:sw_pres(Z(miss(1)-1),lat),lat);
depths = mkgrid(Zext,2*Zext(1)-Zext(2));
T = interp1(Z,T1,depths,'linear');
S = interp1(Z,S1,depths,'linear');

%P = sw_pres(depths,ilat);

N2 = bfrq(S,T,depths,lat);%10^(-6)*ones(32,1);
N2(end) = N2(end-1); %kludge
[Vmode, Hmode, c] = vertmode(N2,depths,n_modes);
%close; close; % close mode plots

% Note: N2, Hmode, Vmode are on Zext = avg1(depths)
% Add a vertically invariant barotropic mode for horizontal velocities
Hmode(:,2:n_modes+1) = Hmode;
Hmode(:,1) = 1/sqrt(size(Hmode,1));
n_modes = n_modes + 1;


%% Do fit

fprintf('\n\n Doing fit.\n');
baroclinic = 1;
barotropic = 0;
taper_param = 0.008;
svd_param = 0.01;

% ignore top because of non wave mode contribution
% NaN's at bottom
d_range = 34:size(U_all,1); % 220m - 3530m
t_range = 1:size(U_all,2);

% mat suffix means 2D'ed grid variable (using repmat)
tmat = repmat(tgrid(t_range),length(d_range),1);
zmat = repmat(z_all(d_range)',1,length(t_range));

% Subsample mode shapes at data depths and combine
zfit = z_all(d_range);
umode = interp1(Zext,Hmode(:,1:end),zfit,'linear');
N2_zfit = interp1(Zext,N2,zfit,'linear')';
% figure; plot(umode',zfit); revz;

ufit = addnan(U_all(d_range,t_range),150);%(detrend(U_all','constant'))';
vfit = addnan(V_all(d_range,t_range),150);%(detrend(U_all','constant'))';
ubfit = nan_detrend(ufit);
vbfit = nan_detrend(vfit);

% Fit all modes (incl barotropic) to data
if barotropic
    [Au_simple Au_taper Au_svd] = fit_1D(ufit,umode,tgrid(t_range),z_all(d_range),taper_param,svd_param,0);
    [Av_simple Av_taper Av_svd] = fit_1D(vfit,umode,tgrid(t_range),z_all(d_range),taper_param,svd_param,0);
end

% Fit baroclinic modes to anomalies
if baroclinic
    [Aub_simple Aub_taper Aub_svd] = fit_1D(ubfit,umode(:,2:n_modes),tgrid(t_range),z_all(d_range),taper_param,svd_param,0);
    [Avb_simple Avb_taper Avb_svd] = fit_1D(vbfit,umode(:,2:n_modes),tgrid(t_range),z_all(d_range),taper_param,svd_param,0);
end

% Reconstruction

if barotropic
    usimple = (umode*Au_simple');
    utaper = (umode*Au_taper');
    usvd = (umode*Au_svd');

    vsimple = (umode*Av_simple');
    vtaper = (umode*Av_taper');
    vsvd = (umode*Av_svd');
end

if baroclinic
    ubsimple = (umode(:,2:end)*Aub_simple');
    ubtaper = (umode(:,2:end)*Aub_taper');
    ubsvd = (umode(:,2:end)*Aub_svd');

    vbsimple = (umode(:,2:end)*Avb_simple');
    vbtaper = (umode(:,2:end)*Avb_taper');
    vbsvd = (umode(:,2:end)*Avb_svd');
end


%% Check mode orthogonality at each time step
check_orth =0;
if check_orth
    
    d_range12 = 34:size(U_all,1);
    
    ufit1 = addnan(U_all(d_range12,t_range),150);%(detrend(U_all','constant'))';
    vfit1 = addnan(V_all(d_range12,t_range),150);%(detrend(U_all','constant'))';

    masku = ~isnan(ufit1);
    maskv = ~isnan(vfit1);
    clear resu resv

    for i=1:length(tgrid)
        % [res] = chk_mode_orth(mode,N2,Z)
        resu(:,:,i) = chk_mode_orth(umode(masku(:,i),:),N2_zfit(masku(:,i)), zfit(masku(:,i)));
        resv(:,:,i) = chk_mode_orth(umode(maskv(:,i),:),N2_zfit(masku(:,i)), zfit(masku(:,i)));
    end
    
    for i=1:length(tgrid)
        imagesc(resu(:,:,i));
        caxis([0 1]);
        colorbar;
        title(i);
        pause(0.01);
    end
end

%% Examine Fit

exam_fit = 1;

if exam_fit
    fprintf('\n\n Examining fit. \n');
    t_range3 = 1:60;

    if barotropic
        plot_amp('u',Au_simple(t_range3,:),Au_taper(t_range3,:),Au_svd(t_range3,:),tgrid(t_range3));
        plot_amp('v',Av_simple(t_range3,:),Av_taper(t_range3,:),Av_svd(t_range3,:),tgrid(t_range3));
    end
    
    if baroclinic
        plot_amp('ub',Aub_simple(t_range3,:),Aub_taper(t_range3,:),Aub_svd(t_range3,:),tgrid(t_range3));
        plot_amp('vb',Avb_simple(t_range3,:),Avb_taper(t_range3,:),Avb_svd(t_range3,:),tgrid(t_range3));
    end

    if barotropic
        plot_energy('u',Au_simple(t_range3,:),Au_taper(t_range3,:),Au_svd(t_range3,:),[0:n_modes-1]);
        plot_energy('v',Av_simple(t_range3,:),Av_taper(t_range3,:),Av_svd(t_range3,:),[0:n_modes-1]);
    end

    if baroclinic
        plot_energy('ub',Aub_simple(t_range3,:),Aub_taper(t_range3,:),Aub_svd(t_range3,:),[1:n_modes-1]);
        plot_energy('vb',Avb_simple(t_range3,:),Avb_taper(t_range3,:),Avb_svd(t_range3,:),[1:n_modes-1]);
    end
end

%% Plot / Analyze Spectra

f_width = 3; % >=2 (3 or 5)
delta_t = 2; % in units of days
d_range2 = 1:10:length(d_range);
t_range2 = 1:60; % time range where full water column is available

% High pass filter everything above 20 day periods.
% " this code was built for making N-day running averages from hourly data, 
% so, if you want to average over 10 timesteps, the second argument should
% be 10./24."
%Au_svd_hi= Au_svd - run_avg2(Au_svd,20./24,1);
%Av_svd_hi= Av_svd - run_avg2(Av_svd,20./24,1);

clear specu specu_taper specu_simple specu_svd specv specv_taper specv_sim ple specv_svd
clear specub specub_taper specub_simple specub_svd specvb specvb_taper specvb_simple specvb_svd

fprintf('\n\n Calculating spectra. \n');
for i=1:length(d_range2)
    %fprintf('\n\n i = %d \n\n', i);
    if barotropic
        [specu(i,:),~, coheramp_out,coherpha_out,freq] = coher(ufit(i,t_range2),ufit(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);    
        [specu_simple(i,:),~, coheramp_out,coherpha_out,freq] = coher(usimple(i,t_range2),usimple(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);    
        [specu_taper(i,:),~, coheramp_out,coherpha_out,freq] = coher(utaper(i,t_range2),utaper(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);
        [specu_svd(i,:),~, coheramp_out,coherpha_out,freq] = coher(usvd(i,t_range2),usvd(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);

        [specv(i,:),~, coheramp_out,coherpha_out,freq] = coher(vfit(i,t_range2),vfit(1,t_range2),delta_t,f_width,'v hi mode1','u hi mode2',0);    
        [specv_simple(i,:),~, coheramp_out,coherpha_out,freq] = coher(vsimple(i,t_range2),vsimple(1,t_range2),delta_t,f_width,'v hi mode1','u hi mode2',0);    
        [specv_taper(i,:),~, coheramp_out,coherpha_out,freq] = coher(vtaper(i,t_range2),vtaper(1,t_range2),delta_t,f_width,'v hi mode1','u hi mode2',0);
        [specv_svd(i,:),~, coheramp_out,coherpha_out,freq] = coher(vsvd(i,t_range2),vsvd(1,t_range2),delta_t,f_width,'v hi mode1','u hi mode2',0);
    end
    
    if baroclinic
        [specub(i,:),~, coheramp_out,coherpha_out,freq] = coher(ubfit(i,t_range2),ubfit(1,t_range2),delta_t,f_width,'u_b hi mode1','u_b hi mode2',0);    
        [specub_simple(i,:),~, coheramp_out,coherpha_out,freq] = coher(ubsimple(i,t_range2),ubsimple(1,t_range2),delta_t,f_width,'u_b hi mode1','u_b hi mode2',0);    
        [specub_taper(i,:),~, coheramp_out,coherpha_out,freq] = coher(ubtaper(i,t_range2),ubtaper(1,t_range2),delta_t,f_width,'u_b hi mode1','u_b hi mode2',0);
        [specub_svd(i,:),~, coheramp_out,coherpha_out,freq] = coher(ubsvd(i,t_range2),ubsvd(1,t_range2),delta_t,f_width,'u_b hi mode1','u_b hi mode2',0);

        [specvb(i,:),~, coheramp_out,coherpha_out,freq] = coher(vbfit(i,t_range2),vbfit(1,t_range2),delta_t,f_width,'u hi mode1','u_b hi mode2',0);    
        [specvb_simple(i,:),~, coheramp_out,coherpha_out,freq] = coher(vbsimple(i,t_range2),vbsimple(1,t_range2),delta_t,f_width,'u_b hi mode1','u_b hi mode2',0);    
        [specvb_taper(i,:),~, coheramp_out,coherpha_out,freq] = coher(vbtaper(i,t_range2),vbtaper(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);
        [specvb_svd(i,:),~, coheramp_out,coherpha_out,freq] = coher(vbsvd(i,t_range2),vbsvd(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);
    end
end

% Average spectra over depth & plot 
if barotropic
    plot_spec('u',specu,specu_simple,specu_taper,specu_svd,freq);
    plot_spec('v', specv,specv_simple,specv_taper,specv_svd,freq);
end

if baroclinic
    plot_spec('u_b',specub,specub_simple,specub_taper,specub_svd,freq);
    plot_spec('v_b', specvb,specvb_simple,specvb_taper,specvb_svd,freq);
end

%% Sensitivity of spectra to fit parameters

sens = 0;
if sens
    taper_param = 0.008;
    svd_param = 0.01;
    str = ['ub ', num2str(svd_param)];
    
    %[Au_simple Au_taper Au_svd] = fit_1D(ufit,umode,tgrid(t_range),z_all(d_range),0.08,0.01,0);
    [Aub_simple Aub_taper Aub_svd] = fit_1D(ubfit,umode(:,2:n_modes),tgrid(t_range),z_all(d_range),taper_param,svd_param,0);
    %usimple = (umode*Au_simple');
    %ubtaper = (umode(:,2:n_modes)*Aub_taper');
    ubsvd = (umode(:,2:n_modes)*Aub_svd');

    for i=1:length(d_range2)
        %fprintf('\n\n i = %d \n\n', i);
        %[specu(i,:),~, coheramp_out,coherpha_out,freq] = coher(ufit(i,t_range2),ufit(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);    
        %[specu_simple(i,:),~, coheramp_out,coherpha_out,freq] = coher(usimple(i,t_range2),usimple(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);    
        %[specub_taper(i,:),~, coheramp_out,coherpha_out,freq] = coher(ubtaper(i,t_range2),ubtaper(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);
        [specub_svd(i,:),~, coheramp_out,coherpha_out,freq] = coher(ubsvd(i,t_range2),ubsvd(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);
    end
    %close
        
    plot_spec(str,specu,[],[],specub_svd,freq);
    plot_energy(str,[],[],Aub_svd,[1:n_modes-1]);
    plot_amp(str,[],[],Aub_svd(t_range2,:),tgrid(t_range2));
end

%% plot code

% ENERGY
% figure;plot([0:n_modes-1],nanmean(Au_simple(t_range2,:).^2),'r'); title('SIMPLE'); set(gca,'XTick',[0:1:n_modes-1]);
% figure;plot([0:n_modes-1],nanmean(Au_taper(t_range2,:).^2),'r'); title('TAPER'); set(gca,'XTick',[0:1:n_modes-1]);
% figure;plot([0:n_modes-1],nanmean(Au_svd(t_range2,:).^2),'r'); title('SVD'); set(gca,'XTick',[0:1:n_modes-1]);

% figure
% levs=[-20:1:20];
% contourf(tgrid,z_all,U_all,levs)
% shading flat
% set(gca,'ydir','reverse')

% Energy plots

% figure;
% hold on
% plot(mean(Au.^2),'r');
% plot(mean(Av.^2),'g');
% plot(mean(At.^2),'b');

% remember first column is barotropic mode
% coher(Au_svd(:,2),Au_svd(:,3),delta_t,f_width,'U mode1',' U mode2',1);
% coher(Avb_svd(:,2),Avb_svd(:,3),delta_t,f_width,'V mode1',' V mode2',1);
% 
% [u1,u2, coheramp_out,coherpha_out,freq] = coher(Au_svd_hi(:,2),Au_svd_hi(:,3),delta_t,f_width,'u hi mode1','u hi mode2',1);
% [v1,v2, coheramp_out,coherpha_out,freq] = coher(Av_svd_hi(:,2),Av_svd_hi(:,3),delta_t,f_width,'V hi mode1','V hi mode2',1);
% 
% figure
% title('Mode 1 u/v ratio');
% semilogy(freq,u1./v1);
% grid on
% 
% figure
% title('Mode 2 u/v ratio');
% semilogy(freq,u2./v2);
% grid on

    %     figure;
    %figure; semilogy(freq,mean(specu)); grid on; title('original u'); ax = axis;
%     semilogy(freq,mean(specub)); 
%     hold on;
%     semilogy(freq,mean(specub_svd),'r'); 
%     grid on; 
%     title(['ub ', num2str(svd_param), ' red - taper']);

    %figure;plot([0:n_modes-1],nanmean(Au_taper(t_range2,:).^2),'r'); title('TAPER'); set(gca,'XTick',[0:1:n_modes-1]);

%% MISC


% d_mmp = 25:500; % range of depths to consider
% d_adcp = 51:116;
% t_range = 1:60;

% Uadcp = ugrid_adcp(d_adcp,t_range);
% Ummp =  ugrid(d_mmp,t_range);
% Vadcp = vgrid_adcp(d_adcp,t_range);
% Vmmp =  vgrid(d_mmp,t_range);

% Ufit = [ugrid_adcp(d_adcp,t_range); ugrid(d_mmp,t_range)];
% Vfit = [vgrid_adcp(d_adcp,t_range); vgrid(d_mmp,t_range)];

% Zadcp = z(d_adcp);
% Zmmp = sw_dpth(pgrid(d_mmp),lat);
% close all;

% t0 = datenum('10-31-2000 05:00:00');
%t1 = tgrid;%t0 + time/24/60;
% plot(mean(At.^2),'b');

%Fw = Vmode(d_range+d,:);
% gap_len = 20; %in units of days
% len = [180]/delta_t; % length of data record (in days) we're searching for
%d = find(Zext == Zin(1))-1;

% Interpolate climatological values to depth grid (for isotherm disp.)
% Tz1 = diff(T)./diff(depths);
% Tz = Tz1(d_range+d);%interp1(Zmid,Tz1,[depths; 800],'linear');
% Ttao = Tgrid(d_range,t_range) - repmat(mean(Tgrid(d_range,t_range),2),1,t_range(end));
% eta = -Ttao./repmat(Tz',t_range(end),1)';

%     Fu_adcp = interp1(Zext,Hmode,Zadcp,'linear');
%     Fu_mmp = interp1(Zext,Hmode,Zmmp,'linear');

%% SVD methods (eta)

% [U1,lam1,V1] = svd(Fw);
% K1 = find(diag(lam1)>0.1.*max(max(lam1)));
% UK1 = U1(:,K1);
% VK1 = V1(:,K1);
% lamK1 = lam1(K1,K1);
% 
% Fpinv1=VK1*(lamK1\UK1');
% At = (Fpinv1*eta)';
