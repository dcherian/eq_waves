function [filtered] = conv_band_pass(in,windows)
    
    if size(in,1) == 1; in = in'; transpose_flag = 1; end
    
    filtered = conv(in,ones(1,windows(1))./windows(1),'same') - ...
        conv(in,ones(1,windows(2))./windows(2),'same');

    filtered = filtered(windows(2):end-windows(2)) - ...
        nanmean(filtered(windows(2):end-windows(2)));
    
    if transpose_flag, filtered = filtered'; end
