function new_bodyt = move_body_z_shift(new_body,flagellum,z_shift)
    temp = flagellum(1,:) - new_body(end,:) + [0,0,z_shift];
    new_bodyt = new_body + repmat(temp,size(new_body,1),1);
end