function [V,H,c] = modes(lon,lat)
	% Load and process data
	
	% ncload('woa05_sal.cdf');
	% ncload('woa05_temp.cdf');
	% sal = salinity;
	% temp = temperature;
	% miss = find(sal == -99.9999008178711);
	% sal(miss) = NaN;
	% 
	% miss = find(temp == -99.9999008178711);
	% temp(miss) = NaN;
	
	if exist('temp') ~= 1 & exist('sal') ~=1
        printf('\n Loading WOA 05 data');
        load woa05.mat;
	end
	
	% Locate (lon,lat)
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
	
	ilon = x(1);
	ilat = y(1);
	Y(ilat)
    X(ilon)
	t = temp(:,ilat,ilon);
	s = sal(:,ilat,ilon);
	p = sw_pres(Z,ilat);
	
	% Plot data
	close all;
	%figure;
	%plot(squeeze(s),Z);
	%xlabel('Salinity (PSU)');
	%ylabel('Z (m)');
	%set(gca,'ydir','reverse');
	%figure;
	%plot(squeeze(t),Z);
	%xlabel('Temperature (C)');
	%ylabel('Z (m)');
	%set(gca,'ydir','reverse');
	
	% Obtain buoyancy frequency
	N2 = bfrq(s,t,Z,ilat);
	[N3,vort,p_ave] = sw_bfrq(s,t,p,ilat); % returns N^2
	
	figure
	err = (N2-N3)./N3 * 100;
	plot(err,Z(1:length(Z)-1));
	xlabel('Percentage error in N2');
	ylabel('Z (m)');
	set(gca,'ydir','reverse');
	
	%N2 = 1*10^-6*ones(32,1);
	
	% Obtain vertical modes
	[V, H, c] = vertmode(N2,Z,3);
end