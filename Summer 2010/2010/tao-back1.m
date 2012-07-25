n = 5; % number of modes

if exist('T_20') ~= 1 
    fprintf('\n Loading TAO & WOA data ...');
    load 'tao_woa.mat';
end

close all;

t0 = datenum('10-31-2000 05:00:00');
t1 = t0 + time/24/60;
lz = length(Z);
Zmid = (Z(1:lz-1) + Z(2:lz))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get vertical modes for location (lon,lat) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('temp') ~= 1 & exist('sal') ~=1
    fprintf('\n Loading WOA 05 data');
    load woa05.mat;
end

% Locate (lon,lat) % ASSUMES EITHER INTEGER OR ___.5 N/E
if isinteger(lon)
    x = find(X < lon+1 & X > lon-1);
else
    x = find(X == lon);
end

if lat < 0 
    y = find(Y < lat+1 & Y > lat-1);
else
    y = find(Y > lat+1 & Y < lat-1);
end

ilon = x(1); ilat = y(1);

% Verify location
fprintf('Actual location = (%.1f W,%.1f N)', 360-X(ilon), Y(ilat));

T = temp(:,ilat,ilon);
S = sal(:,ilat,ilon);
P = sw_pres(Z,ilat);

N2 = bfrq(S,T,Z,ilat);
[Vmode, Hmode, c] = vertmode(N2,Z,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate isotherm displacement %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_int = interp1(Z,T,[depth; 600],'linear');
Tz = diff(T_int)./diff([depth; 600]);

% Find locations of valid data records
%len = [120]; % length of data record (in days) we're searching for
%[ind,num,spillover] = find_gap(T_20,len*24*6);

 a = 410947; b = 451266;
%for i=1:length(len)
%    for j=1:1%num(i)
        %a = ind(j,i); b = ind(j,i) + len(i)*24*6;
        %check_gap(T_20(a:b,:),a,len(i));
        Ttao = detrend(T_20(a:b,6:11),'constant');
        eta = Ttao./repmat(Tz(6:11)',b-a+1,1);
        F = interp1(Zmid,Hmode,depth(6:11),'linear');
        
		% Solve AF'=eta
		 A1 = (F\eta')';
		%for ind=1:length(eta)
   %         nn=eta(ind,:);
    %        A(:,ind)=F'/nn;
    % end
		
		%figure
		%plot(t1(a:b),A)
		%datetick('x','mm/dd/yyyy','keeplimits');
		%ylabel('Modal Amplitude');
		%legend(num2str([n:-1:1]'));
        %   end
        %end

%figure;
%contour(t1(a:b),depth(6:11),T_20(a:b,6:11)',[10 15 20 25]);
%set(gca,'ydir','reverse'); 
%datetick('x','mm/dd/yyyy','keeplimits');

%figure
%MM = repmat(depth(6:11)',b-a+1,1);
%eta1 = eta + MM;
%plot(t1(a:b),eta1);
%set(gca,'ydir','reverse');
%datetick('x','mm/dd/yyyy','keeplimits');
