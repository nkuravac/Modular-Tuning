function new_bodyt = move_body(new_body,flagellum)
     temp = flagellum(1,:) - new_body(end,:);
     new_bodyt = new_body + repmat(temp,size(new_body,1),1);
end