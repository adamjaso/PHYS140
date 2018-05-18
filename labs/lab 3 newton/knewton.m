function value=knewton(E,ecc,angSp,t)
value=-(E-ecc*sin(E)-angSp*t)/(1-ecc*cos(E));