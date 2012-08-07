%% Data
    
load 'data\140W\5N.mat'
%load 'data\modes.mat'

if exist('temp','var') ~= 1 && exist('sal','var') ~=1
    fprintf('\n Loading WOA 05 data');
    load 'data\woa05.mat';
end

%% Calculate theoretical modes

lon = 360-140;
lat = 5;
n_modes = 5;

% Locate (lon,lat)
ilon = find_approx(X,lon,1);
ilat = find_approx(Y,lat,1);

% Verify location
fprintf('\n Actual location = (%.1f W,%.1f N)', 360-X(ilon), Y(ilat));

% Interpolate to denser grid.
Zmode = avg1(Z);
T = temp(:,ilat,ilon);
S = sal(:,ilat,ilon);
dtdz = diff(T)./diff(Z);
%P = sw_pres(depths,lat);

N2 = bfrq(S,T,Z,lat);%10^(-6)*ones(32,1);
[Vmode, Hmode, c] = vertmode(N2,Z,n_modes);
close; close; % close mode plots
Tmode = Vmode .* repmat(dtdz,1,size(Vmode,2));

%% Filter

% Running averages 6 day , 12 day and subtract
dhtavg = conv(dht5n',ones(1,6)./6,'same') - conv(dht5n',ones(1,12)/12,'same');

% FIRST AND LAST TWO VALUES SEEM TO BE CRAP
for i=1:11
    tavg(:,i) = conv(t5n(i,:)',ones(1,6)./6,'same') - conv(t5n(i,:)',ones(1,12)/12,'same');
    tavg(:,i) = tavg(:,i) - nanmean(tavg(:,i));
end

dhtavg = dhtavg - nanmean(dhtavg);

%% Plot 'mode' structure

n_mode = 2;
Imode = dhtavg(240:710,:)\tavg(240:710,:);
ind500 = find_approx(Zmode,500,1);
Imode = Imode./max(Imode);
plot(Imode,depth,'b');
hold on
plot(Tmode(:,n_mode)'./max(Tmode(1:ind500,n_mode)),Zmode,'r');
legend('Inferred mode', 'Theoretical mode');
plot(zeros(size(depth)),depth,'k-')
ylim([depth(1) depth(end)]);
revz