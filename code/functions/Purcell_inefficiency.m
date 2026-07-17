function ineff = Purcell_inefficiency(theory_drag_coeff,bigdir,new_body,freq,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase)
    [~,~,~,avg_sumTf,~,~,~,~,~,~,~,avg_U_net_calc,~]=simulate_bacterium(bigdir,new_body,freq,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
    %theory_drag_coeff is the expected F/U (6*pi*eta*r for a sphere)
    ineff=abs(avg_sumTf(3)*(freq*2*pi)/(theory_drag_coeff*avg_U_net_calc(3).^2));
end