function R = quatern2rotMat(q)
%QUATERN2ROTMAT Converts a quaternion orientation to a rotation matrix
%
%   R = quatern2rotMat(q)
%
%   Converts a quaternion orientation to a rotation matrix.
%
%
%	Date          Author          
%	2021/11/10    Deng zhengxiong    

    R(1,1,:) = q(:,1).^2+q(:,2).^2-q(:,3).^2-q(:,4).^2;
    R(1,2,:) = 2.*(q(:,2).*q(:,3)-q(:,1).*q(:,4));
    R(1,3,:) = 2.*(q(:,1).*q(:,3)+q(:,2).*q(:,4));
    R(2,1,:) = 2.*(q(:,1).*q(:,4)+q(:,2).*q(:,3));
    R(2,2,:) = q(:,1).^2-q(:,2).^2+q(:,3).^2-q(:,4).^2;
    R(2,3,:) = 2.*(q(:,3).*q(:,4)-q(:,1).*q(:,2));
    R(3,1,:) = 2.*(q(:,2).*q(:,4)-q(:,1).*q(:,3));
    R(3,2,:) = 2.*(q(:,1).*q(:,2)+q(:,3).*q(:,4));
    R(3,3,:) = q(:,1).^2-q(:,2).^2-q(:,3).^2+q(:,4).^2;
end

