function [angle,axis] = get_included_angle_from_unit_vector(vec_a, vec_b)
%
%
%	Date          Author         
%	2021/11/11    Deng zhengxiong    

angle_vec = cross(vec_a, vec_b);
angle_sin = norm(angle_vec);
angle_cosin = dot(vec_a, vec_b);
if angle_sin == 0 
    if angle_cosin > 0 % angle is zero
        angle = 0;
        axis = [0,0,0];
        return;
    else
        % angle is 180 degree
        rvec = [0,vec_a(3),-vec_a(2)];
        rvec = rvec / norm(rvec);
        angle = pi;
        axis = rvec;
        return;
    end
else
    if angle_sin > 1
        angle = pi/2;
    else
        angle = asin(angle_sin);
    end
    
    if angle_cosin < 0
        angle = pi - angle;
    end
    
    axis = angle_vec / norm(angle_vec);
end

end
