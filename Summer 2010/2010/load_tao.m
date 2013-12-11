function [] = load_tao(name)
	s = sprintf('data\\cdf\\%s.cdf',char(name));
    fprintf('\n Loading...');
    ncload(s);
	clear('-regexp','QT_*','ST_*');
	m1 = find(T_20 > 100);
	T_20(m1) = NaN;
    cd('C:\dc');
	s1 = sprintf('data\\mat\\%s.mat', name);
    fprintf('\n Saving...');
    save(s1)
    
    a = input('\n Clear workspace (y/n Default y)?; ');
    if isempty(a) | char(toupper(a)) == 'Y' 
        clear;
    end
    