close all
% id 328141809
% 1D ----------------------------------

% values
n_0 = 1.909;
d=9.675*10.^-4;
theta_B = 62.352;
c = 3*10.^8;  % meter / s


lambda_center=1.309*10.^-6;
lambda_values = linspace(lambda_center- 2 * 10.^-9,lambda_center+ 2 * 10.^-9, 1000);  

nu = @(lambda) c./lambda;
n_glass = @(lambda, a) (n_0 - a./(nu(lambda)).^2);
theta_glass = @(lambda, a) asind(sind(theta_B)./n_glass(lambda, a));  % transfer angle
R = @(lambda, a) (sind(theta_B-theta_glass(lambda, a))./(sind(theta_B+theta_glass(lambda, a)))).^2;
delta = @(lambda, a) ((4*180*d*cosd(theta_B)./lambda) .* n_glass(lambda, a));
u = @(lambda, a) (sind(delta(lambda, a)./2)).^2;
transfer = @(lambda, a) ((1-R(lambda, a)).^2)./((1-R(lambda, a)).^2+4*R(lambda, a).*u(lambda, a));

a_values = [10.^26, 10.^27, 10.^28];
nu_FSR_values = zeros(size(a_values));
lambda_FSR_values = zeros(size(a_values));

hold on;
for i = 1:length(a_values)
   transfer_values = transfer(lambda_values, a_values(i));
   plot(lambda_values, transfer_values, 'LineWidth', 2);

   lambda_values_for_FSR = linspace(lambda_center - 100 * 10.^-9,lambda_center + 100 * 10.^-9, 1000);  

   [transfer_peak, lambda_for_max] = findpeaks(transfer_values);
   frquency_values = nu(lambda_values_for_FSR(lambda_for_max));
   nu_FSR_values(i) = mean(abs(diff(frquency_values)));
   lambda_FSR_values(i) = mean(abs(diff(lambda_values_for_FSR(lambda_for_max))));
end

% axes
xlabel('lambda [m]')
ylabel('transfer of intensity')
title('Roy Turgeman - transfer of intensity as a function of the wave length for multiple values of a');

%Display the grid
grid on;

% Legend
legend(cellstr(num2str(a_values', 'a = %-d')));

hold off;

% Display FSR values
disp('nu FSR values for each value of a:');
disp(nu_FSR_values);
disp('lambda FSR values for each value of a:');
disp(lambda_FSR_values);




