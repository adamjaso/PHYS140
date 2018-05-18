%gcdist.m
%This program calculates the total of minimum great circle distances from
%latitude and longitude coordinates of different cities on the surface of
%an (approximately) spherical earth. Angles West and South are to be
%entered as negative values

radiusEarth=6371.01;

%Input
cities=input('Enter the number of cities to be visited: ');
units=input('Enter the units of your coordinates ([d]egrees or [r]adians): ','s');
count=1;
coords=input('Enter first city [lat,long]: ');
while count < cities,
    temp=input('Enter next city [lat,long]: ');
    coords=[coords;temp];
    count=count+1;
end

%Convert to radians
if units == 'd',
    coords = pi / 180 * coords;
end

%Calculation

count=2;
sum = 0;
fprintf('\n');
while count <= cities,
    distance = dist(coords(count-1,:),coords(count,:));
    sum = sum + distance;
    fprintf('The distance between city %d and city %d is %5.5f km\n',count-1,count,distance);
    count=count+1;
end
fprintf('The total distance of this trip was %5.5f km\n',sum);