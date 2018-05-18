%This script calculates the area for a radius that user is prompted
radius = input('Please enter the radius:');
% call calcarea to compute the area of the cicle for this radius 
area = calcarea(radius)
fprintf('For a circle with a radius of %.2f,',radius);
fprintf(' the area is %.4f',area);