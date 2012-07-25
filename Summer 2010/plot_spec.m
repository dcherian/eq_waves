function [] = plot_spec(str,var,simple,taper,svd,freq)

    figure; 
    semilogy(freq,mean(var),'k','LineWidth',1.35); 
    title(str);
    hold on;grid on; 
    if ~isempty(simple), semilogy(freq,mean(simple),'r','LineWidth',1.35); end;
    if ~isempty(taper), semilogy(freq,mean(taper),'g','LineWidth',1.35); end;
    if ~isempty(svd), semilogy(freq,mean(svd),'b','LineWidth',1.35); end;
    if isempty(simple) || isempty(taper) || isempty(svd)
        legend('Original','Depends on arguments');
    else
        legend('Original', 'simple', 'taper', 'svd');
    end