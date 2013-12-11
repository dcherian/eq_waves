% Checks whether the vertical modes specified in Mode are orthogonal.
% Returns sum over all depths.
%       [res] = chk_mode_orth(Mode,N2,Z)
%           Mode -> dimensions Depths x Modes
%           N2 -> Buoyancy frequency squared in (s^-2)
%           Z -> depth grid

function [res] = chk_mode_orth(mode,N2,Z)
    if size(mode,1) ~= length(Z)
        if size(mode,2) == length(Z)
            fprintf('\n Flipping Mode. Check');
            mode = mode';
        else
           fprintf('\n Incorrect dimensions for mode shapes. Check \n');
           return
        end
    end
    
    n_modes = size(mode, 2);
    
    for i=1:n_modes
        for j=1:n_modes
            int1 = mode(:,i).*mode(:,j).*N2;
            if size(int1) ~= size(Z)
                Z = Z';
            end
            res(i,j) = sum(avg1(int1).* diff(Z));
        end
        if res(i,i) ~= 0
            res(i,:) = res(i,:)./(res(i,i));
        end
    end    