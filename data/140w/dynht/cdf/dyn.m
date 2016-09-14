files = ls('*.cdf');

for i=1:size(files,1)
    [vars atts dims] = ncdfread(files(i,:));
    dyn(i,:) = squeeze(vars.DYN_13);
end