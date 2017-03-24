function [range] = findCommonTimeRange(timedht, timetemp)

    start = find(timetemp == timedht(1));
    stop = find(timetemp == timedht(end));
    range = start:stop;

    assert(all(timetemp(range) == timedht));
end
