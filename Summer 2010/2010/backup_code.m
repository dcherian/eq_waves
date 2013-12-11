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
