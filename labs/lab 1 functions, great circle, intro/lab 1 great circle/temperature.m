%temperature.m
%This program calculates the temperature of a planet given the latitude, lat, the luminosity, L,
%the albedo, A, and the distance from star to planet, r.
 
albedoes=[.1,.4];
emissivities=[.97,1.0,.8,.5];
lat = input('Enter the latitude: ');
L = 1;
r = 1;
 
for R=1:2,
    for C=1:4,
        A=albedoes(1,R);
        emissivity=emissivities(1,C);
        temp = 395 * ( ( L * ( 1 - A ) * cos( lat ) ) / ( emissivity * r ^ 2 ) )^.25;
        fahr = 9 / 5 * ( temp - 273.15 ) + 32;
        fprintf('The temperature of the planet is %10f Kelvin and %10f Fahrenheit\n',temp,fahr);
        fprintf('     albedo = %5.5f     emissivity = %5.5f     \n',A,emissivity);
    end
end
