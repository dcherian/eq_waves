% Calculates equatorial wave  spectrum from Eriksen (1982) pg. 1220

function [Ek,Ey,Eg] = eriksen_spec(c)
    % k -> Kelvin
    % y -> Yanai / mixed rossby-gravity
    % g -> gravity
    
    for i=1:length(c)
        Ek(i) = 0.3 * 10^9 * ((1+c(6)/c(i))^(-2))./((sum(1+c(6)./c))^(-2));
        Ey(i) = 1.0 * 10^9 * ((1+c(6)/c(i))^(-2))./((sum(1+c(6)./c))^(-2));
        Eg(i) = 1.0 * 10^9 * (c(i).*(1+c(6)/c(i))^(-2))./((sum(c.*(1+c(6)./c)))^(-2));
    end
    
    figure
    semilogy([1:length(c)],Ek,'r','LineWidth',1.5);
    hold on
    semilogy([1:length(c)],Ey,'g','LineWidth',1.5);
    semilogy([1:length(c)],Eg,'b','LineWidth',1.5);
    legend('Kelvin','Yanai','Gravity');
    grid 
    
    figure
    plot([1:length(c)],Ek,'r','LineWidth',1.5);
    hold on
    plot([1:length(c)],Ey,'g','LineWidth',1.5);
    plot([1:length(c)],Eg,'b','LineWidth',1.5);
    legend('Kelvin','Yanai','Gravity');
    grid on
    