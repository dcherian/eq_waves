% loads data from cdf files and pre-processes them
PROCESS_CDF = 1;
cd('C:\dc');

if PROCESS_CDF == 1    
	D = dir('data\cdf\*.cdf');
	m = size(D);
	
	for i=1:size(D)
        load_tao(D.name(i));
	end
else
    D = dir('data\mat\*.mat');
    for i=1:size(D)
        s = sprintf('data\mat\%s',D.name(i));
        fprintf('\n Loading ... %s', D.name(i));
        load s;
    end
    
    if ~exist('temp') | ~exist('sal')
        load 'dat\mat\woa.mat';
    end
end