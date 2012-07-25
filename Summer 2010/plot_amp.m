function [] = plot_amp(str,simple,taper,svd,xax)
    if ~isempty(simple), figure; plot(xax,simple); title(['SIMPLE ', str]); end;
    if ~isempty(taper), figure; plot(xax,taper); title(['TAPER ', str]); end
    if ~isempty(svd), figure; plot(xax,svd); title(['SVD ', str]); end