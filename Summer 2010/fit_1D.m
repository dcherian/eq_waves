function [A_simple, A_taper, A_svd, res] = fit_1D(data,mode,tgrid,zgrid,taper_param,svd_param,plot)

    mask = ~isnan(data);
    n_modes = size(mode,2);

    for i=1:length(tgrid)
        A_simple(i,:) = simple_ls(data(mask(:,i),i),mode(mask(:,i),:));
        A_taper(i,:) = taper_ls(data(mask(:,i),i),mode(mask(:,i),:),taper_param); % last argument = alpha^2 
        A_svd(i,:) = svd_fit(data(mask(:,i),i),mode(mask(:,i),:),svd_param); % last argument factor -> find(diag(lam)>factor.*max(max(lam)));
    end

    % Reconstruct
    usimple = (mode*A_simple').*fillnan(double(mask),0);
    utaper = (mode*A_taper').*fillnan(double(mask),0);
    usvd = (mode*A_svd').*fillnan(double(mask),0);

    % Plots
     if(plot)
%         figure;plot([0:n_modes-1],nanmean(A_simple.^2),'r'); title('SIMPLE'); set(gca,'XTick',[0:1:n_modes-1]);
%         figure;plot([0:n_modes-1],nanmean(A_taper.^2),'r'); title('TAPER'); set(gca,'XTick',[0:1:n_modes-1]);
%         figure;plot([0:n_modes-1],nanmean(A_svd.^2),'r'); title('SVD'); set(gca,'XTick',[0:1:n_modes-1]);
% 
%         figure; plot(tgrid,A_simple); title('SIMPLE');
%         figure; plot(tgrid,A_taper); title('TAPER');
%         figure; plot(tgrid,A_svd); title('SVD');

        figure
        levs=[-20:1:20];
        contourf(tgrid,zgrid,data,levs)
        shading flat; colorbar
        set(gca,'ydir','reverse')
        caxis([-20 20]);

%         figure
%         levs=[-20:1:20];
%         contourf(tgrid,zgrid,usimple,levs);
%         title('SIMPLE');
%         shading flat; colorbar
%         set(gca,'ydir','reverse')
%         caxis([-20 20]);

        figure
        levs=[-20:1:20];
        contourf(tgrid,zgrid,utaper,levs);
        title('TAPER');
        shading flat; colorbar
        set(gca,'ydir','reverse')
        caxis([-20 20]);
% 
%         figure
%         levs=[-20:1:20];
%         contourf(tgrid,zgrid,usvd,levs);
%         title('SVD');
%         shading flat; colorbar
%         set(gca,'ydir','reverse')
%         caxis([-20 20]);
    end