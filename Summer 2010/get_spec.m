function [specinitu, specfitu_simple, specfitu_taper, specfitu_svd, freq ] = get_spec()

    [specinitu(i,:),~, coheramp_out,coherpha_out,freq] = coher(ufit(i,t_range2),ufit(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);    
    [specfitu_simple(i,:),~, coheramp_out,coherpha_out,freq] = coher(usimple(i,t_range2),usimple(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);    
    [specfitu_taper(i,:),~, coheramp_out,coherpha_out,freq] = coher(utaper(i,t_range2),utaper(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);
    [specfitu_svd(i,:),~, coheramp_out,coherpha_out,freq] = coher(usvd(i,t_range2),usvd(1,t_range2),delta_t,f_width,'u hi mode1','u hi mode2',0);