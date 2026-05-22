function ineff = Purcell_inefficiency(bigdir,new_body,freq,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase)
    [~,~,avg_sumFf,avg_sumTf,~,~,~,~,~,~,~,avg_U_net_calc,~]=simulate_bacterium(bigdir,new_body,freq,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
    ineff=abs(avg_sumTf(3)*(freq*2*pi)/(avg_sumFf(3)*avg_U_net_calc(3)));
end