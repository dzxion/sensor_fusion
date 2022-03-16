function pos = map_projection_project( lat, lon, alt )
%   lat : deg
%   lon : deg
%   alt : cm
%	Date          Author         
%	2021/11/16    Deng zhengxiong

global ToRad ToDeg;
rEarth = 637139300; % µØÇò°ë¾¶ cm
lat = lat * ToRad;
lon = lon * ToRad;
lat_rad_distance = lat - lat(1);
lon_rad_distance = lon - lon(1);
x = lat_rad_distance * rEarth;
x = x';
y = lon_rad_distance.*(rEarth * cos(lat)); 
y = y';
z = alt - alt(1);
z = -z';
pos = [x,y,z];
end