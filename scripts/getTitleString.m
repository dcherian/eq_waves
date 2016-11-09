% Given lon, lat return formatted title string
function [str] = getTitleString(lon, lat)

    if lon > 0
        lonstr = 'E';
    else
        lonstr = 'W';
    end

    if lat < 0
        latstr = 'S';
    else
        latstr = 'N';
    end

    str = sprintf('(%3.1f%s, %1.1f%s)', ...
                  abs(lon), lonstr, ...
                  abs(lat), latstr);
end