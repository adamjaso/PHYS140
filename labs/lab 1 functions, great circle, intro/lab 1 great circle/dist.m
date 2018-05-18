function distance = dist(coord1,coord2)
%This function calculates the great circle distance between two sets of
%latitude and longitude coordinates

%theta = latitude
%phi = longitude
radiusOfEarth = 6371.01;
theta1 = coord1(1,2);
theta2 = coord2(1,2);
phi1 = coord1(1,1);
phi2 = coord2(1,1);

ANGLE=acos( sin( theta1 ) * sin( theta2 ) + cos( theta1 ) * cos( theta2 ) * cos( phi1 - phi2 ));

distance = radiusOfEarth * ANGLE;