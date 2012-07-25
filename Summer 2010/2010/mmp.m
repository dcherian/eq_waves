n_modes = 15; % number of modes
n_moor = 1;

if exist('ugrid','var') ~= 1 
     fprintf('\n Loading MMP data ...');
     load 'mmp_atlantic.mat';
end

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
% Zin -> input pressure grid
% Zext -> input grid extended to surface and bottom (bottom as in WOA'05)
% depths -> grid such that midpoints are at Zin

T1 = temp(:,ilat,ilon);
S1 = sal(:,ilat,ilon);

miss = find(isnan(T1) == 1);
if isempty(miss), miss = length(Z); end

Zin = sw_dpth(pgrid,lat)';
dP = pgrid(2)-pgrid(1);
Zext = sw_dpth(dP:dP:sw_pres(Z(miss(1)-1),lat),lat);
depths = mkgrid(Zext,2*Zext(1)-Zext(2));
T = interp1(Z,T1,depths,'linear');
S = interp1(Z,S1,depths,'linear');
%P = sw_pres(depths,ilat);

N2 = bfrq(S,T,depths,lat);%10^(-6)*ones(32,1);
N2(end) = N2(end-1);
[Vmode, Hmode, c] = vertmode(N2,depths,n_modes);
%close; close; % close mode plots

% N2, Hmode, Vmode are on Zext

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate spectrum %%
%%%%%%%%%%%%%%%%%%%%%%%%%

f_width = 3; 
delta_t = 3;%10/60/24; % in units of days
gap_len = 20; %in units of days
len = [180]/delta_t; % length of data record (in days) we're searching for
d_range = 25:500; % range of depths to consider
t_range = 1:60;
d = find(Zext == Zin(1))-1;

% Interpolate climatological values to depth grid (for isotherm disp.)
Tz1 = diff(T)./diff(depths);
Tz = Tz1(d_range+d);%interp1(Zmid,Tz1,[depths; 800],'linear');
Ttao = Tgrid(d_range,t_range) - repmat(mean(Tgrid(d_range,t_range),2),1,t_range(end));
eta = -Ttao./repmat(Tz',t_range(end),1)';

Fu = Hmode(d_range+d,:);%interp1(Zmid,Hmode,depths(d_range),'linear');
Fw = Vmode(d_range+d,:);

% Simple Lesat Squares
% Au = ugrid(d_range,t_range)'/Fu';
% Av = vgrid(d_range,t_range)'/Fu';
% At = eta'/Fw';

% Use SVD method (u & v)
[U,lam,V] = svd(Fu);
K = find(diag(lam)>0.1.*max(max(lam)));
UK = U(:,K);
VK = V(:,K);
lamK = lam(K,K);

Fpinv=VK*(lamK\UK');
Au = (Fpinv*ugrid(d_range,t_range))';
Av = (Fpinv*vgrid(d_range,t_range))';

% SVD methods (eta)
[U1,lam1,V1] = svd(Fw);
K1 = find(diag(lam1)>0.1.*max(max(lam1)));
UK1 = U1(:,K1);
VK1 = V1(:,K1);
lamK1 = lam1(K1,K1);

Fpinv1=VK1*(lamK1\UK1');
At = (Fpinv1*eta)';

% coher(Au(:,1),Au(:,2),delta_t,f_width,'mode1','mode2',1);
% coher(Av(:,11),Av(:,12),delta_t,f_width,'mode1','mode2',1);
% coher(At(:,11),At(:,12),delta_t,f_width,'mode1','mode2',1);

% Energy plots
figure;
hold on
plot(mean(Au.^2),'r');
plot(mean(Av.^2),'g');
plot(mean(At.^2),'b');

% some TAO data are on a diff. Z-grid. so interpolate to standard depth grid
% if isequal(depth,depths) | isequal(depth',depths)
%     Tfill = T_20;
% else
% 	fprintf('\n Interpolating to standard depths.');
% 	Tfill = (interp1(depth,T_20',depths,'linear'))'; 
% end

%Tfill = fill_gap(Tfill(:,d_range),'linear',gap_len/delta_t); % fill in 'gap_len' day gaps
%[ind,num,spillover] = find_gap(Tfill,len);

% a = 410947; b = 451266;
% for i=1:length(len)
%     for j=1:1%num(i)
        %a = ind(j,i); b = ind(j,i) + len(i);
      %  fprintf('\n Double checking data in loop j=%d.', j);
        %check_gap(Tfill(a:b,:),1,len(i));
        
        % Vertical Isotherm displacement
        
        
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
% %         end
%     end % for j     
%     fprintf('\n Averaging spectra.');
%    % spectra_avg(:,:,i) = mean(spectra(:,:,:,i),3);
%  end % for i		

 
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
