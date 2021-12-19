function q = axisAngle2quatern(axis, angle)
%AXISANGLE2QUATERN Converts an axis-angle orientation to a quaternion
%
%   q = axisAngle2quatern(axis, angle)
%
%   Converts and axis-angle orientation to a quaternion where a 3D rotation
%   is described by an angular rotation around axis defined by a vector.
%
%
%	Date          Author          
%	2021/11/11    Deng zhengxiong    

    q0 = cos(angle./2);
    q1 = axis(:,1)*sin(angle./2);
    q2 = axis(:,2)*sin(angle./2);
    q3 = axis(:,3)*sin(angle./2);
    q = [q0 q1 q2 q3];
end

