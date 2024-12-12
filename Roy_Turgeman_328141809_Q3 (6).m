% 328141809

% -----------------------------
% Warm up quetion
close all;
% 1B
R = 0.05;
L = 0.2;
N = 200;

Delta_x = L/N;
Delta_y = L/N;

x_values = -L/2 : Delta_x : L/2 - Delta_x;
y_values = -L/2 : Delta_y : L/2 - Delta_y;

matrix = zeros(length(x_values), length(y_values));

for i = 1:length(x_values)
    for j = 1:length(y_values)
        matrix(i, j) = circ(x_values(i), y_values(j), R);
    end
end

% 1C

figure;
imagesc(x_values, y_values, matrix);
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title(['1C : Roy Turgmean - circ function with R=', num2str(R), ' [m]']);
axis xy;
axis image;

% 1D
fft_matrix = fftshift(fft2(matrix));  % fft calculation

% Frequency range, where f_x,f_y are the frequency values
f_x_values = -1/(2*Delta_x) : 1/L : 1/(2*Delta_x) - 1/L;
f_y_values = -1/(2*Delta_y) : 1/L : 1/(2*Delta_y) - 1/L;

% Absolute value of the fft of circ plot , with imagesc
figure;
imagesc(f_x_values, f_y_values, abs(fft_matrix))
xlabel('f_x [cycles/m]');
ylabel('f_y [cycles/m]');
title(['1D : Roy Turgmean - absolute value of fft of circ using imagesc with R=', num2str(R), ' [m]']);
colorbar;
axis xy;
axis image;

% Absolute value of the fft of circ plot , with surf
figure;
surf(f_x_values, f_y_values, abs(fft_matrix));
camlight left;
lighting phong;
shading interp;
title(['1D : Roy Turgmean - absolute value of fft of circ using surf with R=', num2str(R), ' [m]']);
xlabel('f_x [cycles/m]');
ylabel('f_y [cycles/m]');
zlabel('Absolute value of fft of circ');


%----------------------------------------
% Question 2

% Define the function parameters
R_final = R;

lambda = 1.129*10^-6;
lambda_final = lambda;

z_0 = 5 * R^2/lambda;  % 2A
z_final = z_0;

f_x_values = -1/(2*Delta_x) : 1/L : 1/(2*Delta_x) - 1/L;
f_y_values = -1/(2*Delta_y) : 1/L : 1/(2*Delta_y) - 1/L;

x_values = lambda_final*z_0 * f_x_values;
y_values = lambda_final*z_0 * f_y_values;

[x, y] = meshgrid(x_values, y_values);

z = ((R_final)^4/(lambda_final * z_final)^2) * (jinc((R_final/(lambda_final*z_final) * sqrt(x.^2+y.^2)))).^ 2;  % 2B

% 2C - plot the function z
figure;
imagesc(x_values, y_values, real(z)); % plot only the real part of the function z
title('2B : Roy Turgmean - plot of the intensity at z=z_0');
xlabel('x [m]');
ylabel('y [m]');
colorbar;
axis xy;
axis image;

% 2D - plot the 2D cut:
my_2D_cut = real(z(N/2 + 1, :));
[right_point, left_point, fwhm, x_interpolation, interpolated_cut] = CalcFWHM(my_2D_cut, x_values);

figure;
hold on;
plot(x_interpolation, interpolated_cut);
x_range = linspace(left_point, right_point, 100);
plot(x_range, max(interpolated_cut)/2 * ones(100));
xlabel('y [m]');
ylabel('intesity [W]');
title('2D : Roy Turgeman - 2D Cut of the intensity where x = 0');
disp(['Fwhm in 2D: ', num2str(fwhm)]);
hold off;

% ------------------------------------------------
% Question 3

z_final = z_0/(50);

x_values = -L/2 : Delta_x : L/2-Delta_x;
y_values = -L/2 : Delta_y : L/2-Delta_y;
[x, y] = meshgrid(x_values, y_values);
u1 = circ(x,y,R);

% 3C - z=z_0/50
[u_2, x_prop] = propFresnel(u1, L, lambda, z_0/50);
intensity = abs(u_2).^2;

figure;
imagesc(x_prop, x_prop, intensity);
xlabel('f_x [cycles/m]');
ylabel('f_y [cycles/m]');
colorbar;
title(['Roy Turgmean - Intensity plot in fresnel diffraction where z/z_0 = ', num2str(1/50)]);
axis xy;
axis image;

my_2D_cut = intensity(N/2 + 1, :);
[right_point, left_point, fwhm, x_interpolation, interpolated_cut] = CalcFWHM(my_2D_cut, x_prop);  % find fwhm

% plot the 2D cut
figure;
hold on;
plot(x_interpolation, interpolated_cut);
x_range = linspace(left_point, right_point, 100);
plot(x_range, max(interpolated_cut)/2 * ones(100));
xlabel('y [m]');
ylabel('intesity [W]');
title(['Roy Turgeman - 2D Cut of the intensity in frensel diffraction where x=0 and z/z_0 = ', num2str(1/50)]);
hold off;
disp(['Fwhm in 3C: ', num2str(fwhm)]);


% 3D
%figure;
[u_2, x_prop] = propFresnel(u1, L, lambda, z_0);
intensity = abs(u_2).^2;

figure;
imagesc(x_prop, x_prop, intensity);
xlabel('f_x [cycles/m]');
ylabel('f_y [cycles/m]');
colorbar;
title(['Roy Turgmean - Intensity plot in fresnel diffraction where z/z_0 = ', num2str(1)]);
axis xy;
axis image;

my_2D_cut = intensity(N/2 + 1, :);
[right_point, left_point, fwhm, x_interpolation, interpolated_cut] = CalcFWHM(my_2D_cut, x_prop);  % find fwhm

% plot the 2D cut
figure;
hold on;
plot(x_interpolation, interpolated_cut);
plot(x_interpolation, interpolated_cut);
x_range = linspace(left_point, right_point, 100);
plot(x_range, max(interpolated_cut)/2 * ones(100));
xlabel('y [m]');
ylabel('intesity [W]');
title(['Roy Turgeman - 2D Cut of the intensity in frensel diffraction where x=0 and z/z_0 = ', num2str(1)]);
disp(['Fwhm in 3D: ', num2str(fwhm)]);
hold off;


% 3 E
z_values = [0.01*z_0, 0.1*z_0, 0.5*z_0, z_0, 2*z_0, 10*z_0];
fwhm_array = zeros(1, length(z_values)); 

figure('Position', [500, 600, 1600, 600]); % Adjust position and size as needed
annotation('textbox', [0.5, 0.9, 0.1, 0.1], 'String', 'Roy Turgeman : 3D plot above 2D cut plot', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'LineStyle', 'none');

for i=1:length(z_values)
    curr_z = z_values(i);
    subplot(2, length(z_values), i);

    [u_2, x_prop] = propFresnel(u1, L, lambda, curr_z);
    intensity = abs(u_2).^2;

    imagesc(x_prop, x_prop, intensity);
    xlabel('f_x [cycles/m]');
    ylabel('f_y [cycles/m]');
    colorbar;
    title(['z/z_0 = ', num2str(curr_z/z_0)]);
    %title(['Roy Turgmean - Intensity plot in fresnel diffraction where
    %z/z_0 = ', num2str(curr_z/z_0)]);
    axis xy;
    axis image;

    % find fwhm
    my_2D_cut = intensity(N/2, :);
    [right_point, left_point, fwhm, x_interpolation, interpolated_cut] = CalcFWHM(my_2D_cut, x_prop);
    
    % plot the 2D cut
    subplot(2, length(z_values), i+length(z_values));
    hold on;
    plot(x_interpolation, interpolated_cut);
    plot(x_interpolation, interpolated_cut);
    x_range = linspace(left_point, right_point, 100);
    plot(x_range, max(interpolated_cut)/2 * ones(100));
    xlabel('x [m]');
    ylabel('intesity [(V/m)^2]');
    title(['z/z_0 = ', num2str(curr_z/z_0)]);
    hold off;
    %title(['Roy Turgeman - 2D Cut of the intensity in frensel diffraction where y=0 and z/z_0 = ', num2str(curr_z/z_0)]);

    fwhm_array(i) = fwhm;
end


% plot fwhm as a function of z/z_0
figure;
plot(z_values/z_0, fwhm_array, 'o-', 'LineWidth', 1.5);
xlabel('z/z_0');
ylabel('FWHM (m)');
title('Roy Turgeman - FWHM as a function of z/z_0');
grid on;


function [right_point, left_point, fwhm, x_interpolation, interpolated_cut] = CalcFWHM(cut, x_prop)
    x_interpolation = linspace(min(x_prop), max(x_prop), 100); % Interpolation range
    interpolated_cut = interp1(x_prop, cut, x_interpolation, 'spline');

    [max_cut, max_index] = max(interpolated_cut);
    half_max_cut = max_cut / 2;
    
    % find where the cut first drops below half_max_cut:
    % before the peak and after the peak
    left_index = find(interpolated_cut(1:max_index) < half_max_cut, 1, 'last');
    right_index = find(interpolated_cut(max_index:end) < half_max_cut, 1) + max_index - 1;
    
    % x values at the indices
    left_point = x_interpolation(left_index);
    right_point = x_interpolation(right_index);
    
    % Calculate FWHM
    fwhm = right_point - left_point;
end

