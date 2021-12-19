function euler = quatern2euler(q)
%QUATERN2EULER Converts a quaternion orientation to ZYX Euler angles
%
%   euler = quatern2euler(q)
%
%   Converts a quaternion orientation to ZYX Euler angles where roll is a
%   rotation around X, pitch around Y and yaw around Z.
%
%
%	Date          Author         
%	2021/11/09    Deng zhengxiong    

    roll = atan2( 2.*(q(:,1).*q(:,2)+q(:,3).*q(:,4)) , 1-2.*(q(:,2).^2+q(:,3).^2) );
    pitch = asin( 2.*(q(:,1).*q(:,3)-q(:,2).*q(:,4)) );
    yaw = atan2( 2.*(q(:,1).*q(:,4)+q(:,2).*q(:,3)) , 1-2.*(q(:,3).^2+q(:,4).^2) );

    euler = [roll(:,1) pitch(:,1) yaw(:,1)];
end

