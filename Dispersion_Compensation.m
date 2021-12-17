function CDC = Dispersion_Compensation(r_srrc,sps,sps_new,fs_new,Ddcf,L2)
    r=resample(r_srrc,sps_new,sps);%down to sps=2
    %CDC fiber in dsp area
    f_new=fs_new*(-0.5:1/length(r):0.5-1/length(r));
    CDC=fiber(r,f_new,Ddcf,L2);
end

