% Returns & plots normalized variable & vertical co-ordinate (WKB approximation)
%       [Vn,Zn] = norm_wkb(var1,N2,Z)
%           Vn -> normalized variable 'var1'
%           Zn -> normalized vertical co-ordinate Z
%           N2 -> Buoyancy frequency squared (at each Z)
%           var1 -> variable to normalize. Each column should correspond to
%           Zmid.

function [Vn,Zn] = norm_wkb(var1,N2,Z)
	
	N0 = 3/3600; % 3cph
    N1 = sqrt(N2)/N0;
    
    % If filled in at bottom. Assuming N2 is at midpoints of Z grid
%     if length(N1)-length(Z) == -1
%         N1 = N1(1:end-1);
%     end

    dZn = (N1(1:end-1)+N1(2:end))/2.*diff(Z);
    Zn(2:length(Z)) = cumsum(dZn);
    Zn(1) = Z(1);
    
    Vn = var1./repmat(sqrt(N1),1,size(var1,2));
    figure;
    plot(Vn,Zn);
    revz;