d_range=3:10;
Tfill = (interp1(depth,u_1205',depths,'linear'))';

Tfill = fill_gap(Tfill(:,d_range),'linear',gap_len/delta_t); % fill in 'gap_len' day gaps
[ind,num,spillover] = find_gap(Tfill,len);