n_modes = 20; % number of modes
n_moor = 1;

if exist('T_20','var') ~= 1 
    fprintf('\n Loading TAO & WOA data ...');
    load 'tao_woa.mat';
end
% else
%     % clear all variables except those in tao_woa.mat
% %     clear i j var2 var1;
% %     var2 = char(who);
% %     var1 = char(who('-file','tao_woa.mat'));
% %     
% %     for i=1:length(var2)
% %         for j=1:length(var1)
% %             if strcmp(var2(i,:),var1(j,:))
% %                 break;                
% %             end
% %             clear(var2(i,:));
% %         end
% %     end
% end

%close all;

t0 = datenum('10-31-2000 05:00:00');
t1 = t0 + time/24/60;

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
ilon = x(1); ilat = y(1);

% Verify location
fprintf('\n Actual location = (%.1f W,%.1f N)', 360-X(ilon), Y(ilat));

% Interpolate to denser grid.
Zmid = avg1(Z);
T = temp(:,ilat,ilon);
S = sal(:,ilat,ilon);
%P = sw_pres(depths,lat);

N2 = bfrq(S,T,Z,lat);%10^(-6)*ones(32,1);
[Vmode, Hmode, c] = vertmode(N2,Z,n_modes);
close; close; % close mode plots

% N2, Hmode, Vmode are on Zmid

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate spectrum %%
%%%%%%%%%%%%%%%%%%%%%%%%%

f_width = 5; 
delta_t = 10/60/24; % in units of days
gap_len = 20; %in units of days
len = [100]/delta_t; % length of data record (in days) we're searching for
depths = [1 25 50 75 100 125 150 200 250 300 500]'; % standard depths
d_range = [2:5]; % range of depths to consider

% Interpolate climatological values to standard depth grid (for isotherm
% displacement)
Tz1 = diff(T)./diff(Z);
Tz = interp1(Zmid,Tz1,[depths; 800],'linear');
F = interp1(Zmid,Hmode,depths(d_range),'linear');
 
% some TAO data are on a diff. Z-grid. so interpolate to standard depth grid
if isequal(depth,depths) || isequal(depth',depths)
    Tfill = T_20;
else
	fprintf('\n Interpolating to standard depths.');
	Tfill = (interp1(depth,T_20',depths,'linear'))'; 
end

Tfill = fill_gap(Tfill(:,d_range),'linear',gap_len/delta_t); % fill in 'gap_len' day gaps
[ind,num,spillover] = find_gap(Tfill,len);

% a = 410947; b = 451266;
for i=1:length(len)
    for j=1:1%num(i)
        a = ind(j,i); b = ind(j,i) + len(i);
        fprintf('\n Double checking data in loop j=%d.', j);
        check_gap(Tfill(a:b,:),1,len(i));
        
        % Vertical Isotherm displacement
        Ttao = detrend(Tfill(a:b,:),'constant');
        eta = -Ttao./repmat(Tz(d_range)',len(i)+1,1);
	    F1=F; A = eta/F1';
        
        % Calculate integrated hydrostatic pressure
        %Ttao = Tfill(a:b,:);
        %Stao = 35*ones(len(i)+1,length(d_range));
        %rho = sw_dens(Stao, Ttao, repmat(sw_pres(depths(d_range)',lat),len(i)+1,1));
        %p1 = fliplr(cumsum(fliplr(9.81*(rho(:,1:end-1)+ rho(:,2:end))/2.*repmat(diff(depths(2:11)'),len(1)+1,1)),2));
        %F1 = [F(1:9,:); zeros(1,n_modes)];
        %p1(:,10) = 0;
        %A = p1/F1';
        
        % Solve AF'=eta
        %fprintf('\n Carrying out least squares fit.');
        
        coher(A(:,1,j),A(:,2,j),delta_t,f_width,'mode1','mode2',1);
%       coher(mode_amp(1,:,j),mode_amp(2,:,j),delta_t,f_width,'mode1','mode2',1);

%         for k=1:2:n_modes
%             if(k == n_modes) % only for odd n_modes
%                 mode1 = sprintf('mode %d', n_modes-1);
%                 mode2 = sprintf('mode %d', n_modes);
%                 [spectra(:,n_modes-1,j,i),spectra(:,n_modes,j,i), coheramp_out,coherpha_out,freq] = coher(A(:,n_modes-1,j),A(:,n_modes,j),delta_t,f_width,mode1,mode2,0);
%                 continue    
%             end
%             mode1 = sprintf('mode %d', k);
%             mode2 = sprintf('mode %d', k+1);
%             [spectra(:,k,j,i),spectra(:,k+1,j,i), coheramp_out,coherpha_out,freq] = coher(A(:,k,j),A(:,k+1,j),delta_t,f_width,mode1,mode2,0);
%         end
    end % for j     
    fprintf('\n Averaging spectra.');
   % spectra_avg(:,:,i) = mean(spectra(:,:,:,i),3);
 end % for i		

 
 %%%%%%%%%%%%%%%%%%%
 %%%% MISC CODE %%%%
 %%%%%%%%%%%%%%%%%%%

% Plot eta at depths
%disp_plot(eta,t1(a:b),depths(d_range));
%datetick('x','mm/dd/yyyy','keeplimits');

% r= randn(b-a+1,10)*10;
% eta = eta+r;

%       eta=eta';
% 		Finv=inv(F'*F)*F';
% 		for n=1:length(eta)
%          etai=eta(:,n);
%          ahati=Finv*etai;
%          mode_amp(:,n)=ahati;
% 		end

% Plot modal amplitude
% figure
% plot(t1(a:b),A(:,1,1));
% hold on;
% %plot(t1(a:b),A(:,1,1))
% datetick('x','mm/dd/yyyy','keeplimits');
% ylabel('Modal Amplitude');
% legend(num2str([1:n]'));
        
%for k=1:length(eta)
%    nn=eta(k,:);
%    A(k,:,j)=nn/F';
%end
         
%[spectry_out,spectrz_out, coheramp_out,coherpha_out,freq] = coher(A(:,1,1),A(:,2,1),10/60/24,5,'mode 1','mode 2');

% Plot temperature contours
% figure;
% contour(t1,depths,T_20',[10 15 20 25 30]);
% set(gca,'ydir','reverse'); 
% set(gca,'XAxisLocation','top');
% datetick('x','mm/dd/yyyy','keeplimits');
