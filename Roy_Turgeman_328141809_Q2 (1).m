close all;
% id 328141809

% 2B ----------------------
R_1 = 41 * 10^-3;
R_2 = 32 * 10^-3;

lambda_values = 0.4:0.01:1;  % the range of lambda
s = @(lambda) lambda.^2;

n_BK7 = @(lambda) sqrt(1 + (((1.03961212 * s(lambda)) ./ (s(lambda) - 0.00600069867)) ...
    + ((0.231792344 * s(lambda)) ./ (s(lambda) - 0.0200179144)) ...
    + ((1.01046945 * s(lambda)) ./ (s(lambda) - 103.560653))));

%disp(n(0.6934))  % needs to be 1.5132

f = @(lambda) R_1 * R_2 / (R_1 + R_2) * 1 ./ (n_BK7(lambda) - 1);
f_values = f(lambda_values);

figure;
plot(lambda_values, f_values, 'LineWidth', 2)  % plot
xlabel('lambda [\mum]');  
ylabel('f [m]');
title('Roy Turgeman - graph of the focus as a function of the wavelength in micro meter');
grid on;


% 2C2 ------------------

% values
L_1 = 70*10^-3;
R = 39*10^-3;
lambda_1 = 0.559; % micro meter

n_lens = n_BK7(lambda_1);
disp(n_lens);

% ABCD matrices
first_free_space_matrix = @(x) [1, x; 0, 1];
lens_matrix = [1, 0; 2*(1-n_lens)/R, 1];
second_free_space_matrix = @(x) [1, x-L_1; 0, 1];

% compute L_2
syms x
total_ABCD_matrix = second_free_space_matrix(L_1+x) * lens_matrix * first_free_space_matrix(L_1);
L_2 = solve(total_ABCD_matrix(1,2) == 0, x);
L_2 = double(L_2);
disp(L_2);

% compute the vector which defines the ray

values_in_lin_space = 1000;
x_case1 = linspace(0, L_1, values_in_lin_space);
x_case2 = linspace(L_1, L_1+L_2, values_in_lin_space);

vec_case1 = zeros(2, length(x_case1)); 
vec_case2 = zeros(2, length(x_case2));

theta_values = -pi/5 : pi/25 : pi/5 ;

figure;
hold on;

for theta = theta_values
    input_vector = [1; theta];

    % Calculate vec for x_case1
    for i = 1:length(x_case1)
        current_matrix = first_free_space_matrix(x_case1(i));
        vec_case1(:, i) = current_matrix * input_vector;
    end

    % lens - before
    last_vec_case1 = vec_case1(:, length(x_case1));
    vec_after_lens = lens_matrix * last_vec_case1;

    % Calculate vec for x_case2
    for i = 1:length(x_case2)
        current_matrix = second_free_space_matrix(x_case2(i));
        vec_case2(:, i) = current_matrix * vec_after_lens;
    end

    combined_vector = [vec_case1, vec_case2];

    % Plot the first coordinate of the vector
    x_positions_combined = [x_case1, x_case2];
    plot(x_positions_combined, combined_vector(1, :), 'LineWidth', 2);
end

hold off;

xlabel('x (m)');
ylabel('r (m)');
title('Plot of the height of the ray (r) for Multiple Theta Values - Roy Turgeman');
legend(arrayfun(@(theta) sprintf('Theta = %.2f', theta), theta_values, 'UniformOutput', false));
xline(L_1, 'r', 'LineWidth', 2, 'Label', 'L1');
xline(L_1 + L_2, 'b', 'LineWidth', 2, 'Label', 'L1 + L2'); 
grid on;


% Q2C3 ------------------------

% values
L_1 = 70*10^-3;
R = 39*10^-3;
lambda_1 = 0.559; % micro meter
lambda_2 = 0.545; % micro meter

% ABCD matrices
first_free_space_matrix = @(x) [1, x; 0, 1];
lens_matrix = @(lambda) [1, 0; 2*(1-n_BK7(lambda))/R, 1];
second_free_space_matrix = @(x) [1, x-L_1; 0, 1];

% compute L_2
syms x
total_ABCD_matrix = second_free_space_matrix(L_1+x) * lens_matrix(lambda_1) * first_free_space_matrix(L_1);
L_2 = solve(total_ABCD_matrix(1,2) == 0, x);
L_2 = double(L_2);

% compute the vector which defines the ray

values_in_lin_space = 1000;
x_case1 = linspace(0, L_1, values_in_lin_space);
x_case2 = linspace(L_1, L_1+L_2, values_in_lin_space);

vec_case1 = zeros(2, length(x_case1)); 
vec_case2 = zeros(2, length(x_case2));

h_values = -1 :0.1 : 1;

figure;
hold on;

lambda_array_values = [lambda_1, lambda_2];
for i = 1:length(lambda_array_values)
    curr_lambda = lambda_array_values(i);

    for h = h_values
        input_vector = [h; 0];
    
        % Calculate vec for x_case1
        for j = 1:length(x_case1)
            current_matrix = first_free_space_matrix(x_case1(j));
            vec_case1(:, j) = current_matrix * input_vector;
        end
    
        % lens - before
        last_vec_case1 = vec_case1(:, length(x_case1));
        vec_after_lens = lens_matrix(curr_lambda) * last_vec_case1;

    
        % Calculate vec for x_case2
        for k = 1:length(x_case2)
            current_matrix = second_free_space_matrix(x_case2(k));
            vec_case2(:, k) = current_matrix * vec_after_lens;
        end
    
        combined_vector = [vec_case1, vec_case2];
    
        % Plot the first coordinate of the vector
        x_positions_combined = [x_case1, x_case2];
        color = 'r';
        if i == 2
            color = 'b';
        end
        plot(x_positions_combined, combined_vector(1, :), 'LineWidth', 2, 'Color', color);
        hold on;
    end
end

hold off;

xlabel('x (m)');
ylabel('r (m)');
title('Plot of the height or the ray (r) for Multiple h Values (red lambda_1, blue lambda_2) - Roy Turgeman');
legend(arrayfun(@(h) sprintf('h = %.2f', h), h_values, 'UniformOutput', false));
xline(L_1, 'r', 'LineWidth', 2, 'Label', 'L1');
xline(L_1 + L_2, 'b', 'LineWidth', 2, 'Label', 'L1 + L2'); 
grid on;

% Set custom y-axis ticks
yticks(-1:0.1:1);


% Q2D -------------------------------
% Q2DB
R = 39*10.^-3; % mm
lambda_1 = 0.409;  % micro meter
lambda_2 = 0.644;  % micro meter
d_1 = 1.1 * 10.^-3;
d_2 = 2.76 * 10.^-3;

lambda_values = 0.4:0.01:1;  % micro meter

s = @(lambda) lambda.^2;

% first graph
n_F2 = @(lambda) sqrt(1 + (((1.34533359 * s(lambda)) ./ (s(lambda) - 0.00997743871)) ...
    + ((0.209073176 * s(lambda)) ./ (s(lambda) - 0.0470450767)) ...
    + ((0.937357162 * s(lambda)) ./ (s(lambda) - 111.886764))));

% find R_3 - Q2DA

R_1 = R;
R_2 = -R;

syms x;
free_space = @(d) [1, d ; 0, 1];
tau_R_1 = @(lambda) [1,0 ; (1-n_BK7(lambda))./((n_BK7(lambda).*R_1)), 1./n_BK7(lambda)];
tau_R_2 = @(lambda) [1,0 ; (n_BK7(lambda)-n_F2(lambda))./((n_F2(lambda).*R_2)), n_BK7(lambda)./n_F2(lambda)];
tau_R_3 = @(lambda) [1,0; (n_F2(lambda)-1)./x, n_F2(lambda)];

Total_ABCD_Matrix = @(lambda) tau_R_3(lambda) * free_space(d_2) * tau_R_2(lambda) * free_space(d_1) * tau_R_1(lambda);  % not piece wise multiplication
First_Matrix = Total_ABCD_Matrix(lambda_1);
Second_Matrix = Total_ABCD_Matrix(lambda_2);

R_3 = solve(First_Matrix(2,1) == Second_Matrix(2,1), x);
R_3 = double(R_3);
disp(R_3);

% Q2DB

new_tau_R_3 = @(lambda) [1,0; (n_F2(lambda)-1)./R_3, n_F2(lambda)];
new_Total_ABCD_Matrix = @(lambda) new_tau_R_3(lambda) * free_space(d_2) * tau_R_2(lambda) * free_space(d_1) * tau_R_1(lambda);  % not piece wise multiplication

f_values = zeros(0, length(lambda_values));
for i=1:length(lambda_values)
    curr_lambda = lambda_values(i);
    curr_total_matrix = new_Total_ABCD_Matrix(curr_lambda);
    f_values(i) = -1./curr_total_matrix(2,1);
end

figure;
hold on;
plot(lambda_values, f_values, 'LineWidth', 2)  % plot

lambda_min = 0.4;
lambda_max = 1;

% find avg values
f_avg = mean(f_values);
n_BK7_avg = mean(n_BK7(lambda_values));
effective_R = 2 * f_avg * (n_BK7_avg - 1);

new_f_values = zeros(0, length(lambda_values));

new_tau_R_1_avg_case = @(lambda) [1,0 ; (1-n_BK7(lambda))./((n_BK7(lambda).*effective_R)), 1./n_BK7(lambda)];
new_tau_R_2_avg_case = @(lambda) [1,0 ; (n_BK7(lambda)-1)./(-effective_R), n_BK7(lambda)];
total_matrix_avg_case = @(lambda)  new_tau_R_2_avg_case(lambda) * free_space(d_1) * new_tau_R_1_avg_case(lambda);

for i = 1:length(lambda_values)
    curr_lambda = lambda_values(i);
    curr_total_matrix = total_matrix_avg_case(curr_lambda);
    new_f_values(i) = -1./curr_total_matrix(2,1);

end

plot(lambda_values, new_f_values, 'LineWidth', 2)  % plot

% Labeling the axes
xlabel('lambda [\mum]');
ylabel('f [m]');
title('Roy Turgeman - graphs of the focus as a function of the wavelength in micro meter');

% Display the grid
legend('show');
grid on;
hold off;


% Q2DC
syms h1 h2;
new_tau_R_3 = @(lambda) [1,0; (n_F2(lambda)-1)./R_3, n_F2(lambda)];

new_total_ABCD_Matrix = @(lambda) new_tau_R_3(lambda) * free_space(d_2) * tau_R_2(lambda) * free_space(d_1) * tau_R_1(lambda);
last_total_ABCD_Matrix = @(lambda) free_space(h2) * new_total_ABCD_Matrix(lambda) * free_space(h1);

h1_array = zeros(1, length(lambda_values));
h2_array = zeros(1, length(lambda_values));
index = 1;

for lambda=lambda_values
    syms h1;
    syms h2;

    last_total_ABCD_Matrix = @(lambda) free_space(h2) * new_total_ABCD_Matrix(lambda) * free_space(h1);
    Curr_last_ABCD_Matrix = last_total_ABCD_Matrix(lambda);
    
    first_eq = Curr_last_ABCD_Matrix(1,1) == 1;
    second_eq = Curr_last_ABCD_Matrix(2,2) == 1;
    third_eq = Curr_last_ABCD_Matrix(1,2) == 0;
    
    solutions = solve([first_eq, third_eq], [h1, h2]);
    % Checks: play with 1,2 and 2,3 and 1,3 -> same solutions!
    h1_array(index) = double(solutions.h1);
    h2_array(index) = double(solutions.h2);

    index = index + 1;
end

figure;

hold on;
plot(lambda_values, h1_array, 'LineWidth', 2);
plot(lambda_values, h2_array, 'LineWidth', 2);
xlabel('lambda [\mum]');
ylabel('h_1 , h_2 [m]');
title('Roy Turgeman - graphs of the principal planes as a function of the wavelength in micro meter');

grid on;
legend('h1', 'h2'); 
hold off;

final_h1 = mean(h1_array);
final_h2 = mean(h2_array);

disp(final_h1);
disp(final_h2);
