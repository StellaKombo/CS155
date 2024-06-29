clear all;close all;clc;
% Define the trajectory
% time of the simulation
dt = 0.01;
tsim_end = 100; % sampling time in secs
tsim = 0:dt:tsim_end;

t_trajectory = tsim;

% Define the obstacles position
x_obstacle1 = 0.99816;
y_obstacle1 = 3.873;
r_obstacle1 = 0.5;
r_ego = 0.2;

r_total1 = r_obstacle1 + r_ego;

%% Step 2 - Define the trajectory
t_traj = tsim;
figure8_size = 3.0;
figure8_period = 20;

x_traj_ref = figure8_size*cos(2*pi*t_traj/figure8_period);
y_traj_ref = 0.5*figure8_size*sin(4*pi*t_traj/figure8_period);

dot_x_traj_ref = -2*pi*figure8_size*sin(2*pi*t_traj/figure8_period)/figure8_period;
dot_y_traj_ref = 0.5*4*pi*figure8_size*cos(4*pi*t_traj/figure8_period)/figure8_period;

ddot_x_traj_ref = -(2*pi)^2*figure8_size*cos(2*pi*t_traj/figure8_period)/(figure8_period^2);
ddot_y_traj_ref = -0.5*(4*pi)^2*figure8_size*sin(4*pi*t_traj/figure8_period)/(figure8_period^2);

for i = 1:length(t_traj)
    theta_traj_ref(i) = atan2(dot_y_traj_ref(i),dot_x_traj_ref(i));
end
for i = 2:length(t_traj)
    theta_traj_ref(i) = fun_unwrap(theta_traj_ref(i-1),theta_traj_ref(i));
end
for i = 1:length(t_traj)
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

    curr_x = [x_traj_ref(1); y_traj_ref(1); theta_traj_ref(1)];

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
    Q_cost = diag([10,10,1]); 
    R_cost = diag([1,1]);

    % Optimizer and solver
    for k = 2:n_horizon  
        objective = objective + (x_var{k} - r_var{k})' * Q_cost * (x_var{k} - r_var{k});
    end

    for k = 1:(n_horizon - 1)
        objective = objective + u_var{k}' * R_cost * u_var{k};
        constraints = [constraints, x_var{k + 1} == Ad_val{k} * x_var{k} + Bd_val{k} * u_var{k}];
        constraints = [constraints, A_ineq1{k}' * x_var{k+1}(1:2) <= b_ineq1{k}];
    end
 
    % Define the parameters going into the controller
    parameters_in = {x_var{1}, r_var{:}, Ad_val{:}, Bd_val{:}, A_ineq1{:}, b_ineq1{:},};
    solutions_out = {[x_var{:}],[u_var{:}]};
    controller_straps = optimizer(constraints, objective,sdpsettings('solver','gurobi','gurobi.qcpdual',1, 'verbose',1,'debug',1),parameters_in,solutions_out);

    counter = 1;

    % Estimation of the non-linear model
    curr_x_bar = [x_traj_ref(1:n_horizon);y_traj_ref(1:n_horizon);theta_traj_ref(1:n_horizon)];
    u_bar = [u_1_traj_ref(1:n_horizon - 1); u_2_traj_ref(1:n_horizon - 1)];

    initial_position_trajectory_mpc = zeros(length(tsim), 4); % Assuming 4D position for ancilliary controller
    % Create a figure outside the loop for combined trajectory
    h_combined = figure;

while (counter < (length(tsim) - 1 - n_horizon))
        input_list = {};
        input_list{1} = curr_x; % initialize the robot position to be the first position in the reference trajectory        sub_counter = 2; 
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
        
        % Define the matrix of the obstacles
        A_obstacle_1 = []; % (xbar - xobs)'*(xbar - xobs)
        b_obstacle_1 = []; % robs* ||xbar - xobs||

        for l = 2:n_horizon
            position_bar = curr_x_bar(1:2,l - 1);                         
            A_obstacle_1(:,l - 1) = (position_bar - [x_obstacle1;y_obstacle1])' * (position_bar - [x_obstacle1;y_obstacle1]);
            b_obstacle_1(l - 1) =  r_total1 * norm(position_bar - [x_obstacle1;y_obstacle1],2);
        end

        for i = 1: n_horizon -1
            input_list{sub_counter} = A_obstacle_1(:,i);
            sub_counter = sub_counter + 1;
        end

        for i = 1: n_horizon -1
            input_list{sub_counter} = -b_obstacle_1(i);
            sub_counter = sub_counter + 1;
        end

        [solutions, diagnostics] = controller_straps{input_list};

        if diagnostics == 1
            error('The problem is infeasible');
        end 
        x_qp_solution = solutions{1};
        u_qp_solution = solutions{2}; % determine the u_mpc from the optimization
        u_mpc = u_qp_solution;
        u_current_point = u_qp_solution(:,1); % u at the current trajectory
        x_mpc = x_qp_solution; % position obtained when the mpc kicks in
        
        % % Define the ancilliary controller lqr that kicks in after the mpc
        % % to make sure within the horizon there is some reduction of the
        % % oscillations
        Ts = 0.001; % define the time
        [Kd, ~, ~] = lqrd(A_approx, B_approx, Q_cost, R_cost,Ts); 
        difference_traj = initial_position - x_mpc;
        u_lqr = Kd * difference_traj;
        u_total = u_current_point + u_lqr;  % this is now the closed loop controller

        %% Step 6 - Forward simulation with the combined new controller
        options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
        [t, qt_ode] = ode45(@(t, qt) actual_dynamics_ugv_ct(t, qt, u_total),[tsim(counter),tsim(counter + 1)],initial_position);
        contoller_path_output = qt_ode;
        initial_position = qt_ode(end,:)'; % the last position of the new total controller is the new position
      
        [t_qp, xt_qp] = ode45(@(t, x) dynamics(x, clfQP(x)), [0, 20], x0_nlp);
        initial_position_trajectory_(counter + 1, :) = initial_position;
     
        initial_position_trajectory_qp_track(counter + 1, :) = xt_qp(1,:);
        x0_nlp = initial_position;
 
        %% Plotting stuff - Update the controller trajectory plot data
        figure(h_combined);
   
        % Plotting the fixed trajectory
        plot(x_traj_ref, y_traj_ref, '--', 'LineWidth', 2, 'DisplayName', 'Reference Trajectory');
        hold on;

        % Plot the trajectroy followed by the controller based on the MPC
        % optimization
        plot(initial_position_trajectory_mpc(:, 1), initial_position_trajectory_mpc(:, 2), 'r*', 'LineWidth', 2, 'DisplayName', 'MPC Controller Trajectory');
        hold on;

        plot(initial_position_trajectory_qp_track(:, 1), initial_position_trajectory_qp_track(:, 2), 'go', 'LineWidth', 2, 'DisplayName', 'NL-QP Controller Trajectory');

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

        counter = counter + 1;
end


