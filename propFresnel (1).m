function [u_2, x_prop] = propFresnel(u1, L, lambda, z)
% propagation – according to Fresnel
% assumes same x and y side lengths and
% uniform sampling
% u1 - source plane field
% L - source plane side length (FOV)
% lambda - wavelength
% z - propagation distance
% u2 - observation plane field
% x_prop – axis in observation plane field

i=1i;
k=2*pi/lambda;
N=size(u1,1);

Delta_x = L/N;
Delta_y= L/N;

x_values = -L/2 : Delta_x : L/2-Delta_x;
y_values = -L/2 : Delta_y : L/2-Delta_y;
[x, y] = meshgrid(x_values, y_values);

coeff = (exp(i*k*z)/(lambda*z*i)) * exp((i*k/(2*z))*(x.^2+y.^2));
field=u1 .* exp((i*k/(2*z)) .* (x.^2+y.^2));
transform = F(field);

displayed_field = coeff .* transform;
f_x_prop = -1/(2*Delta_x) : 1/L : 1/(2*Delta_x) - 1/L;
x_prop = f_x_prop * (lambda*z);
u_2 = displayed_field;

end







