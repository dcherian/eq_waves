%% make combined plot
% makes huge plot with a subplot for each point on the TAO array and plots
% theoretical and inferred mode 

load Imode.mat

hfig = figure; clf;
set(hfig,'Position',[0 0 1600 900]);
hash = githash([mfilename('fullpath') '.m']);
insertAnnotation(['plot_modes: ' hash]);

xlimits = [-0.2 1.2];
ylimits = [0 500];
fontSize = [20 24 28];

nlon = length(modes.lon);
nlat = length(modes.lat);

hax = packfig(nlon,nlat);

for mm=1:nlon
    for nn=1:nlat 
        ind500 = find_approx(modes.zTmode,500,1);
        if isempty(modes.Imode{mm,nn}), continue; end
        subplot_index = sub2ind([nlon nlat],mm,nn);

        axes(hax(subplot_index));
        hold on

        % temp std
        plot(data.tstd{mm,nn}./max(data.tstd{mm,nn}),data.depth{mm,nn},'k');

        % 0 mean flow mode
        plot(squeeze(abs(modes.Tmode(mm,nn,:))) ...
             ./ max(abs(modes.Tmode(mm,nn,1:ind500))), ...
             modes.zTmode,'r'); %'Color',0.6*[1 1 1]); %

        % inferred mode from TAO data
        plot(modes.Imode{mm,nn},modes.depth{mm,nn}, 'bo-', ...
             'LineWidth', 1);

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

% remove locations without data
hax(4).delete;

hax(46).YAxis.Visible = 'off';
hax(46).XLabel.String = '140W';
axes(hax(46)); beautify(fontSize);
hax(46).XAxis.Color = [1 1 1];
hax(46).XLabel.Color = hax(45).XLabel.Color;

%%

figure(hfig)
[ax,h] = suplabel(['Red = Theoretical Mode, Blue = Inferred Mode, ' ...
                   'Black = Std. of temp at depth'],'t');
ax.YLabel.String = 'Z (m)';
set(h,'FontSize',fontSize(3));
%export_fig('-nocrop','-r150','../images/combined.png');
