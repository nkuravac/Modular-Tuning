function Omega_vecs_tmp = angular_vel_in_direction(U,direction,r_vec,sflag)

Omega_vecs_tmp = zeros(size(U));
for iv = 1:size(U,1)
    %iv
    if iv == sflag
        Omega_vecs_tmp(iv,:) = [0 0 0];
    else
        Omega_vec_dir = cross(direction,r_vec(iv,:));
        tmp = norm(Omega_vec_dir);
        if tmp < 1e-10
            Omega_vecs_tmp(iv,:) = [0 0 0];
        else
            Omega_vecs_tmp(iv,:) = dot(U(iv,:),Omega_vec_dir)/(tmp^2)*direction;
        end
    end
end