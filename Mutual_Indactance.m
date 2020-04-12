clc;

c = 3*10^8;
R2 = 2;

syms phi x y z I r
f = r*I/c * sin(phi) / sqrt(z^2 + (r*cos(phi) - x)^2 + (r*sin(phi) - y)^2);
A = int(f,phi,[0 2*pi]);
disp(A);