clear all;close all;clc;
% Define the trajectory
% time of the simulation
dt = 0.02;

% Define the obstacles position
x_obstacle1 = 2.0;
y_obstacle1 = 1.5;
r_obstacle1 = 0.5;


% Define an initial buffer
r_s = 0.2; % robot radius
r_delta_t = 0.01;
r_ego = r_s + r_delta_t;

r_total1 = r_obstacle1 + r_ego;


%% Step 2 - Define the trajectory
t_trajectory = 0:dt:100;
figure8_size = 3.0;
figure8_period = 20;

% Define each point in the trajcetory with their respective time value

x_traj_ref = figure8_size * cos(2*pi*t_trajectory/figure8_period);
y_traj_ref = 0.5 * figure8_size * sin(4*pi*t_trajectory/figure8_period);

dot_x_traj_ref = -2*pi*figure8_size*sin(2*pi*t_trajectory/figure8_period)/figure8_period;
dot_y_traj_ref = 0.5*4*pi*figure8_size*cos(4*pi*t_trajectory/figure8_period)/figure8_period;

ddot_x_traj_ref = -(2*pi)^2*figure8_size*cos(2*pi*t_trajectory/figure8_period)/(figure8_period^2);
ddot_y_traj_ref = -0.5*(4*pi)^2*figure8_size*sin(4*pi*t_trajectory/figure8_period)/(figure8_period^2);


for i = 1:length(t_trajectory)
    theta_traj_ref(i) = atan2(dot_y_traj_ref(i),dot_x_traj_ref(i));
end
for i = 2:length(t_trajectory)
    theta_traj_ref(i) = fun_unwrap(theta_traj_ref(i-1),theta_traj_ref(i));
end
for i = 1:length(t_trajectory)
    if (sin(theta_traj_ref(i)) == 0)
        u_1_traj_ref(i) = dot_x_traj_ref(i)/cos(theta_traj_ref(i));
    else
        u_1_traj_ref(i) = dot_y_traj_ref(i)/sin(theta_traj_ref(i));
    end
    if u_1_traj_ref(i)~=0
        u_2_traj_ref(i) = (ddot_x_traj_ref(i)*sin(theta_traj_ref(i)) - ddot_y_traj_ref(i)*cos(theta_traj_ref(i)))/(-u_1_traj_ref(i));
    else
       u_2_traj_ref(i) = 1.0/(1.0 + theta_traj_ref(i)^2); 
    end

end

%% Step 3 - MPC
% define the mpc parameters
number_states = 3;
number_controls = 2;
n_horizon = 20;

% Define YAMLIPS decision variables: squared real valued matrices of u and x
u_var = sdpvar(repmat(number_controls, 1, n_horizon - 1), repmat(1, 1, n_horizon - 1));
x_var = sdpvar(repmat(number_states, 1, n_horizon), repmat(1, 1, n_horizon));
Ad_val = sdpvar(repmat(number_states, 1, n_horizon - 1), repmat(number_states, 1, n_horizon-1));
Bd_val = sdpvar(repmat(number_states, 1, n_horizon - 1), repmat(number_controls, 1, n_horizon-1));

% define the reference variable
r_var =  sdpvar(repmat(number_states,1,n_horizon),repmat(1, 1, n_horizon));

% Define the inequalities for the obstacles
A_ineq1 = sdpvar(repmat(2, 1, n_horizon-1), repmat(1,1, n_horizon-1));
b_ineq1 = sdpvar(repmat(1, 1, n_horizon-1), repmat(1,1, n_horizon-1));

% Objective functions
objective = 0; 
constraints = []; 
Q_cost = diag([10,10,10]); 
% Q_cost = diag([100,100,10]); 
%R_cost = diag([1,1]); % gets stuck near the obstacle
%R_cost = diag([0.1,0.1]); % avoids the obstacle but is able to go through most of the trajectory.Better than the 1,1
R_cost = diag([0.8,0.8]);


% Optimizer and solver
for k = 2:n_horizon  
    objective = objective + (x_var{k} - r_var{k})' * Q_cost * (x_var{k} - r_var{k});
end

for k = 1:(n_horizon - 1)
    objective = objective + u_var{k}' * R_cost * u_var{k};
end

for k = 1:(n_horizon - 1)
    constraints = [constraints, x_var{k + 1} == Ad_val{k} * x_var{k} + Bd_val{k} * u_var{k}]; % make it a variable of time the reference trajectory
    constraints = [constraints, A_ineq1{k}' * x_var{k+1}(1:2) <= b_ineq1{k}];
end

% Define the parameters going into the controller
parameters_in = {x_var{1}, r_var{:}, Ad_val{:}, Bd_val{:}, A_ineq1{:}, b_ineq1{:},};
solutions_out = {[x_var{:}],[u_var{:}]};
controller_straps = optimizer(constraints, objective,sdpsettings('solver','gurobi','gurobi.qcpdual',1, 'verbose',1,'debug',1),parameters_in,solutions_out);

counter = 1;
curr_x = [x_traj_ref(1); y_traj_ref(1); theta_traj_ref(1)];

% Estimation of the non-linear model
curr_x_bar = [x_traj_ref(1:n_horizon);y_traj_ref(1:n_horizon);theta_traj_ref(1:n_horizon)];
u_bar = [u_1_traj_ref(1:n_horizon - 1); u_2_traj_ref(1:n_horizon - 1)];
x_polar_previous = fun_inertial2polar(curr_x, curr_x,curr_x);
pos_history_hist = curr_x;
controller_hist = [];
pos_history_polar_hist = [];


% Define parameters for CP
non_conformity_scores = zeros(1, length(t_trajectory) - 1);
N_score = 100;
epsilon_values = zeros(1, N_score);
initial_epsilon = 0.1; % Initial value for epsilon
learning_rate = 0.01; % Learning rate for adaptive epsilon updating
epsilon_desired = 0.05;
h_combined = figure;

% Define the values of the auxilliary controller
k1 = 0.01;
k2 = 0.001;
k3 = 0.001;

% Implementation of the MPC and PID
while (counter < (length(t_trajectory) - 1 - n_horizon))
        input_list = {};
        input_list{1} = curr_x; % initialize the robot position to be the first position in the reference trajectory
        sub_counter = 2;    
        for i = 1:n_horizon
            input_list{sub_counter} = [x_traj_ref(counter + i); y_traj_ref(counter + i); theta_traj_ref(counter + i)];
            sub_counter = sub_counter + 1;
        end

        for i = 1: n_horizon - 1
            % Linearization around some reference trajectory -- Eq7 Daniel et
            % al. 
            % A_approx = df/dx|xo = [0 0 -sin(x3) * u1; 0 0 cos(x3)*u1; 0 0 0]
            A_approx = [0, 0, -sin(curr_x_bar(3, i)) * u_bar(1, i);0, 0, cos(curr_x_bar(3, i)) * u_bar(1, i); 0, 0, 0];
            % b_approx = f(x)
            % We dont have a C matrix bc this is a simple example
            B_approx = [cos(curr_x_bar(3, i)), 0;sin(curr_x_bar(3,i)) ,0 ;0, 1];
            
            % Discretization using expm
            new_mat = [[A_approx, B_approx];[zeros(number_controls, (number_states + number_controls))]];
            new_mat_exp = expm(new_mat * dt);
            input_list{sub_counter} = new_mat_exp(1:number_states, 1:number_states);
            sub_counter = sub_counter + 1;
        end

        for i = 1: n_horizon - 1
            % Linearization around some reference trajectory -- Eq7 Daniel et
            % al. 
            % A_approx = df/dx|xo = [0 0 -sin(x3) * u1; 0 0 cos(x3)*u1; 0 0 0]
            A_approx = [0, 0, -sin(curr_x_bar(3, i)) * u_bar(1, i);0, 0, cos(curr_x_bar(3, i)) * u_bar(1, i); 0, 0, 0];
            % b_approx = f(x)
            % We dont have a C matrix bc this is a simple example
            B_approx = [cos(curr_x_bar(3, i)), 0;sin(curr_x_bar(3,i)) ,0 ;0, 1];
            
            % Discretization using expm
            new_mat = [[A_approx, B_approx];[zeros(number_controls, (number_states + number_controls))]];
            new_mat_exp = expm(new_mat * dt);
            input_list{sub_counter} = new_mat_exp(1:number_states, number_states + 1 :end);
            sub_counter = sub_counter + 1;
        end
        
        % Define the matrix of the obstacles for SCP
        A_obstacle_1 = [];
        b_obstacle_1 = []; 
       
        for l = 2:n_horizon
            position_bar = curr_x_bar(1:2,l - 1);                         
            A_obstacle_1(:,l - 1) = -(position_bar - [x_obstacle1;y_obstacle1]);
            b_obstacle_1(l - 1) = [x_obstacle1;y_obstacle1]'*(position_bar - [x_obstacle1;y_obstacle1]) + r_total1 * norm(position_bar - [x_obstacle1;y_obstacle1], 2);
        end


        for i = 1: n_horizon - 1
            input_list{sub_counter} = A_obstacle_1(:,i);
            sub_counter = sub_counter + 1;
        end

        for i = 1: n_horizon -1
            input_list{sub_counter} = -b_obstacle_1(i);
            sub_counter = sub_counter + 1;
        end

        [solutions, diagnostics] = controller_straps{input_list};
        x_mpc = solutions{1};  % position obtained when the mpc kicks in
        u_qp_solution = solutions{2}; % determine the u_mpc from the optimization

        x_mpc_now = [x_mpc(1,2); x_mpc(2,2); x_mpc(3,2)];  % first position from the mpc formulation
        x_polar_desired = fun_inertial2polar(curr_x, x_mpc_now, x_polar_previous);
        x_polar_previous = x_polar_desired;
        
        % Define the control inputs with the total controller -- mpc & PID
        v_val = u_qp_solution(1,1) + k1 * x_polar_desired(1) * cos(x_polar_desired(2));
        if (x_polar_desired(2) ~= 0)
           w_val = u_qp_solution(2,1) + k2 * x_polar_desired(2) + k1 * (sin(x_polar_desired(2)) * cos(x_polar_desired(2))/x_polar_desired(2))*(x_polar_desired(2) + k3 * x_polar_desired(3));
        else
            w_val = u_qp_solution(2,1);
        end

        u_total = [v_val;w_val];
    %% Step 6 - Forward simulation with the combined new controller
       [t_ode45, qt_ode45] = ode45(@(t, qt) actual_dynamics_ugv_ct(t, qt, u_total),[t_trajectory(counter), t_trajectory(counter+1)], curr_x);
        contoller_path_output = qt_ode45;
        curr_x = qt_ode45(end,:)';
        u_mpc(:,counter) = u_qp_solution(:,1);
        x_desired_mpc(:,counter) = x_mpc(:,2);
    
        pos_history_polar_hist = [pos_history_polar_hist, x_polar_desired'];
        pos_history_hist = [pos_history_hist, curr_x];
        controller_hist = [controller_hist, u_total];
    
        initial_position_trajectory(counter + 1, :) = curr_x;

        %% Adaptive Conformal Prediction
       % Assuming you have conformity scores for each cell, you can populate the cell array
       % if statement sbout the threshold of the cells are computed then we
       % can ignore.
       
       % Compute non-conformity score for the current time ste -  sliding
       % scale

        current_non_conformity = norm(contoller_path_output(end,:) - contoller_path_output(1,:),2);
        non_conformity_scores(1:end-1) = non_conformity_scores(2:end);
        non_conformity_scores(end) = current_non_conformity;

       % Print non-conformity scores for analysis or tracking
       fprintf('Non-conformity scores: %s\n', mat2str(non_conformity_scores));

       sort_conf_score = sort(non_conformity_scores);
       qcount = ceil((N_score + 1) * (1 - initial_epsilon));
       Z = sort_conf_score(qcount);
       if null(B_approx) == 2
           disp("doesn't exist") % doesnt exist because the matrix is full column rank
       else
            Z_perp = null(B_approx); 
       end
       % Update tube radius based on Z
       % r0 = sqrt(norm(Z, 2)^2 / 4); need to calculate after i find
       % the value of alpha
       r_delta_t = sqrt(norm(Z, 2)^2 / 4); 

       r_ego = r_s + r_delta_t;

       r_total1 = r_obstacle1 + r_ego;
       
       % Apply threshold to get boolean condition
       boolean_condition = current_non_conformity >= Z;
       % Map boolean condition to 0 or 1
       R_binary = double(boolean_condition);
       initial_epsilon = initial_epsilon + learning_rate * (epsilon_desired - R_binary);
       
       % Rewrite the array of the epsilon values
       epsilon_values(1:end-1) = epsilon_values(2:end);
       epsilon_values(end) = initial_epsilon;
         
       fprintf('desired epsilon value: %s\n', mat2str(initial_epsilon));

       fprintf('desired epsilon value: %s\n', mat2str(epsilon_desired));

        % Define the approximations of the u and x variables
        u_bar = [u_qp_solution(:,2:end),u_qp_solution(:,end)];
        curr_x_bar = [x_mpc(:,2:end), x_mpc(:,end) + [cos(x_mpc(3,end)) * u_qp_solution(1,end)*dt; sin(x_mpc(3,end)) * u_qp_solution(1,end) * dt; u_qp_solution(2,end) * dt]]; 

        %% Step 7 - Visualizations
        figure(h_combined);
        % Plotting the fixed trajectory
        plot(x_traj_ref, y_traj_ref, '--', 'LineWidth', 2, 'DisplayName', 'Reference Trajectory');
        hold on;

        % Plot the trajectroy followed by the controller based on the MPC
        % optimization
        plot(initial_position_trajectory(:, 1), initial_position_trajectory(:, 2), 'go', 'LineWidth', 2, 'DisplayName', 'Controller Trajectory');
        
        % Plot the obstacles
        scatter(x_obstacle1, y_obstacle1, 'k', 'LineWidth', 2, 'DisplayName', 'Obstacle 1');
        X1_example =[x_obstacle1, y_obstacle1];
        viscircles(X1_example, r_total1, 'LineStyle', '--', 'Color', 'r', LineWidth=2.0); % Plot safety radius circle
   

        title('Combined Controller Trajectory vs Reference Trajectory');
        xlabel('X - axis');
        ylabel('Y - axis');
        legend('Location', 'Best');
        grid on;
        axis equal;

        % Capture frame
        frames_combined(counter) = getframe(h_combined);
        % Introduce a delay to observe the animation
        pause(0.01);
        % Draw the animation
        drawnow;
        hold off;
        % disp(r_val1);
        % disp(r_s);
        % disp(r_ego);
        % 
        % figure(h_combined2);
        % 
        % % % Plot the adaptive epsilon values after the loop
        % theta = 0:pi/50:2*pi;
        % x_obstacle1_plot = r_obstacle1 * cos(theta) + x_obstacle1;
        % y_obstacle1_plot = r_obstacle1 * sin(theta) + y_obstacle1;
        % x_obstacle2_plot = r_obstacle2 * cos(theta) + x_obstacle2;
        % y_obstacle2_plot = r_obstacle2 * sin(theta) + y_obstacle2;
        % 
        % axis([-5,5 -3 3])
        % for i = 1:5:length(initial_position)
        %             xego_plot = r_ego * cos(theta) + initial_position(1,i);
        %             yego_plot = r_ego * sin(theta) + initial_position(2,i);
        %             plot(x_traj_ref,y_traj_ref,'b');
        %             hold on;
        %             plot(xego_plot,yego_plot,'k-');
        %             hold on;    
        %             plot(x_obstacle1_plot,y_obstacle1_plot,'r-');
        %             hold on;
        %             plot(x_obstacle2_plot,y_obstacle2_plot,'r-');
        %             % Capture frame
        %             frames_combined(counter) = getframe(h_combined2);
        %             pause(0.01);
        %             drawnow;
        %             hold off;
        % 
        % end
        disp(counter)

        counter = counter + 1;
   
end

figure();
% Plot the adaptive epsilon values after the loop
figure;
plot(1:N_score * (length(tsim) - 1 - n_horizon), epsilon_values, 'o-', 'LineWidth', 2);
title('Adaptive Epsilon Evolution');
xlabel('Iteration');
ylabel('Epsilon Value');
grid on;