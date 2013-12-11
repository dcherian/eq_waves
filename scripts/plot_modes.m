%% make combined plot
% makes huge plot with a subplot for each point on the TAO array and plots
% theoretical and inferred mode 

load Imode.mat

figure(2); clf;
set(gcf,'Position',[0 0 1600 900]);

xlimits = [-0.2 1.2];
ylimits = [0 500];
fontSize = [11 13 14];

nlon = length(modes.lon);
nlat = length(modes.lat);

for mm=1:nlon
    for nn=1:nlat 
        ind500 = find_approx(modes.zTmode,500,1);
        if isempty(modes.Imode{mm,nn}), continue; end
        subplot_index = sub2ind([nlon nlat],mm,nn);
        figure(2);
        subplot(nlat,nlon,subplot_index);
        hold on
        % temp std
        plot(data.tstd{mm,nn}./max(data.tstd{mm,nn}),data.depth{mm,nn},'k');% 'Color',0.5*[1 1 1]);
        % 0 mean flow mode
        plot(squeeze(abs(modes.Tmode(mm,nn,:)))./max(abs(modes.Tmode(mm,nn,1:ind500))), ...
            modes.zTmode,'r')%'Color',0.6*[1 1 1]); % 
        % inferred mode from TAO data
        plot(modes.Imode{mm,nn},modes.depth{mm,nn},'bo-');
        linex(0);
        xlim(xlimits); ylim(ylimits);
        set(gca,'YTick',[0 200 400]);
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

%%

% fix missing lon labels
missing_lon = 140;
figure(2)
for ii = 1:length(missing_lon)
    subplot(nlat,nlon,sub2ind([nlon nlat],find(modes.lon == missing_lon(ii)),nlat))
    xlabel('140 W');
    set(gca,'YTickLabel',[]);
    beautify(fontSize);
    set(gca,'YColor',[1 1 1]);
    box off;
end

figure(2)
[ax,h] = suplabel('Red = Theoretical Mode, Blue = Inferred Mode, Black = Std. of temp at depth','t');
set(h,'FontSize',fontSize(3));
export_fig('-nocrop','-r150','../images/combined.png');
