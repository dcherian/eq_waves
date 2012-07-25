% MMP + adcp

n_modes = 15; % number of BAROCLINIC modes
n_moor = 1;

if exist('ugrid_adcp','var') ~= 1 
     fprintf('\n Loading MMP + ADCP data ...');
     load 'mmpa.mat';
end

%% Prepare data
d_mmp = 25:500; % range of depths to consider
d_adcp = 51:116;
t_range = 1:60;

% Uadcp = ugrid_adcp(d_adcp,t_range);
% Ummp =  ugrid(d_mmp,t_range);
% Vadcp = vgrid_adcp(d_adcp,t_range);
% Vmmp =  vgrid(d_mmp,t_range);

Ufit = [ugrid_adcp(d_adcp,t_range); ugrid(d_mmp,t_range)];
Vfit = [vgrid_adcp(d_adcp,t_range); vgrid(d_mmp,t_range)];

Zadcp = z(d_adcp);
Zmmp = sw_dpth(pgrid(d_mmp),lat);
%close all;

%t0 = datenum('10-31-2000 05:00:00');
t1 = tgrid;%t0 + time/24/60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get vertical modes for location (lon,lat) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
close; close; % close mode plots

% Note: N2, Hmode, Vmode are on Zext = avg1(depths)
% Add a vertically invariant barotropic mode for horizontal velocities
Hmode(:,2:n_modes+1) = Hmode;
Hmode(:,1) = 1/sqrt(size(Hmode,1));
n_modes = n_modes + 1;

% Subsample mode shapes at data depths and combine
Fu_adcp = interp1(Zext,Hmode,Zadcp,'linear');
Fu_mmp = interp1(Zext,Hmode,Zmmp,'linear');
Fu = [Fu_adcp; Fu_mmp];
Zfit = [Zadcp, Zmmp];

figure; plot(Fu',Zfit); revz;
% Same as above -> use to check
% figure;
% plot(Fu_adcp',Zadcp);
% revz;
% hold on;
% plot(Fu_mmp', Zmmp);

% Use SVD method (u & v)
[UU,lam,VV] = svd(Fu);
% [UU,lam,VV] = svd(Fu(:,2:end));
K = find(diag(lam)>0.1.*max(max(lam)));
UK = UU(:,K);
VK = VV(:,K);
lamK = lam(K,K); 

Fpinv=VK*(lamK\UK');
% Fpinv = pinv(Fu);
% Fit barotropic mode to full U/V
% Aub = Fpinv*Ufit;
% Avb = Fpinv*Vfit;
% figure; % Energy plot
% hold on
% plot(mean(Aub'.^2),'r');
% plot(mean(Avb'.^2),'g');

% Fit baroclinic modes to anomalies
% Au = Fpinv*detrend(Ufit,'constant');
% Av = Fpinv*detrend(Vfit,'constant');
% figure; % Energy plot
% hold on
% plot(mean(Au'.^2),'r');
% plot(mean(Av'.^2),'g');

%% Plot Spectra

f_width = 3; 
delta_t = 3;%10/60/24; % in units of days

% remember first column is barotropic mode
%coher(Au(:,1),Au(:,2),delta_t,f_width,'mode1','mode2',1);
coher(Aub(:,2),Aub(:,3),delta_t,f_width,'BC mode1',' BC mode2',1);
%coher(At(:,11),At(:,12),delta_t,f_width,'mode1','mode2',1);

figure
imagesc(Ufit);
title('Ufit');
figure
imagesc(Fu*Aub);
title('U reconstructed');

figure
imagesc(Vfit);
title('Vfit');
figure
imagesc(Fu*Avb);
title('V reconstructed');
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

% Simple Least Squares
% Au = ugrid(d_range,t_range)'/Fu';
% Av = vgrid(d_range,t_range)'/Fu';
% At = eta'/Fw';



% SVD methods (eta)
% [U1,lam1,V1] = svd(Fw);
% K1 = find(diag(lam1)>0.1.*max(max(lam1)));
% UK1 = U1(:,K1);
% VK1 = V1(:,K1);
% lamK1 = lam1(K1,K1);
% 
% Fpinv1=VK1*(lamK1\UK1');
% At = (Fpinv1*eta)';

% Energy plots
%figure;
% hold on
% plot(mean(Au.^2),'r');
%plot(mean(Av.^2),'g');
% plot(mean(At.^2),'b');