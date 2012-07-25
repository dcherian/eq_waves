function [] = plot_energy(str,simple,taper,svd,xax)

    if ~isempty(simple), figure;plot(xax,nanmean(simple.^2),'r'); title(['SIMPLE ', str]); set(gca,'XTick',xax); end;
    if ~isempty(taper), figure;plot(xax,nanmean(taper.^2),'r'); title(['TAPER ', str]); set(gca,'XTick',xax); end;
    if ~isempty(svd), figure;plot(xax,nanmean(svd.^2),'r'); title(['SVD ', str]); set(gca,'XTick',xax); end;