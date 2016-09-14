files = ls('*.cdf');

for i=1:size(files,1)
    [vars atts dims] = ncdfread(files(i,:));
    tmp(i,:,:) = squeeze(vars.T_20);
end