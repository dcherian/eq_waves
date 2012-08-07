%% make combined plot

load Imode.mat

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

export_fig('-nocrop','-r150','../images/combined.png');