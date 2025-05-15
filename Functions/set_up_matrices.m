function [Hess,Hess_inv,G_eq,g_eq,h_ineq_w,H_ineq,C_A,A,B,C,x0,nx,nu,ny] = set_up_matrices(model,Q_epsilon,n)
    [A,B,C,x0] = start_stript(model);
    nx = size(A,2);
    ny = size(C,1);
    nu = size(B,2);

    if model == "model1"
        y_max = 10*ones(nx,1); y_min = -y_max;
        u_max = 2 *ones(nu,1); u_min = -u_max;
    elseif model == "model2"
        y_max = 30*ones(nx,1); y_min = -y_max;
        u_max = 10*ones(nu,1); u_min = -u_max;
    else
        y_max = 30*ones(nx,1); y_min = -y_max;
        u_max = 8 *ones(nu,1); u_min = -u_max;
    end
    
    R_epsilon       = Q_epsilon;
    [Hess,Hess_inv] = make_Hessian(nx,Q_epsilon,R_epsilon,n,nu);
     C_A            = make_C_A(A,C,n);
    [H_ineq,G_eq]   = make_H_eq_G_ineq(ny,nu,n,A,B,C);
    [h_ineq_w,g_eq] = make_h_eq_q_ineq(G_eq,y_max,y_min,u_max,u_min,n,nx,nu);
end



function [Hess,Hess_inv] = make_Hessian(nx,Q_epsilon,R_epsilon,n,nu)
    state_horizon = 2*nx*n;
    u_horizon = 2*nu*n;
    
    %Quadratic weight 
    Q_hat = eye(state_horizon) * Q_epsilon;
    R_hat = eye(u_horizon) * R_epsilon;
    Hess = blkdiag(Q_hat,R_hat); %[Q_hat 0; 0 R_hat] matrix

    Q_hat_inv = eye(state_horizon) * 1/Q_epsilon;
    R_hat_inv = eye(u_horizon) * 1/R_epsilon;
    Hess_inv = blkdiag(Q_hat_inv,R_hat_inv); 
end

function C_A = make_C_A(A,C, N)
    %Make the [A;A^2;A^3;...;A^N] vector
    [rows, cols] = size(A);
    A_exp_vector = zeros(rows * N, cols); %
    current_power = A;
    for i = 1:N
        start_idx = (i-1)*rows + 1;
        end_idx = i*rows;
        A_exp_vector(start_idx:end_idx,:) = current_power;
        current_power = current_power*A; 
    end

    %Make the C_hat matrix where C is along the diagonal
    I = eye(N);
    C_diag_matrix = kron(I, C);

    %C_A = C_hat*A_exp_vector
    C_A = C_diag_matrix*A_exp_vector;
end
function [H_ineq,G_eq] = make_H_eq_G_ineq(ny, nu, n, A, B, C)
    % Dimensions
    len_A = size(A, 1);

    % Construct Ia matrix
    v = ones(n - 1, 1); % Vector for subdiagonal
    Ia = eye(n * len_A) - kron(diag(v, -1), A);

    % Construct block diagonal matrices for B_hat and C_hat
    B_diag = kron(eye(n), B);
    C_diag = kron(eye(n), C); 

    % Compute C_b (C_diag * Ia^-1 * B_diag)
    C_b = C_diag / Ia * B_diag;

    % Preallocate G_top
    eye_ny = eye(ny * n); % Precompute for reuse
    G_eq = [eye_ny, -eye_ny, -C_b, C_b];

    % Construct G_ineq
    H_under_x = blkdiag([eye_ny; -eye_ny], [eye_ny; -eye_ny]);
    H_under_u = blkdiag([eye(nu*n); -eye(nu*n)], [eye(nu*n); -eye(nu*n)]);

    H_ineq = blkdiag(H_under_x, H_under_u);
end


function [h_ineq,g_eq] = make_h_eq_q_ineq(G_eq,ymax,ymin,umax,umin,n,nx,nu)
    g_eq          = zeros(size(G_eq,1),1);
    
    y_max = kron(ones(n,1),ymax);
    y_min = kron(ones(n,1),ymin);
    u_max = kron(ones(n,1),umax);
    u_min = kron(ones(n,1),umin);
    x_zero = zeros(nx*n,1);
    u_zero = zeros(nu*n,1);
    h_ineq = [y_max;x_zero;-y_min;x_zero;u_max;u_zero;-u_min;u_zero];
end

function [a,b,C,x0] = start_stript(model)
    folder = "Test_models";
    input_folder = fullfile(folder,model);
    % if model == "model1"
    %     C = eye(2);
    % elseif model == "model2"
    %     %C = [0 0 -0.098 0.269;0 0 0.080 0.327];
    %     C = eye(4);
    % elseif model == "model3"
    %     %C = [0 0 -0.098 0.269;0 0 0.080 0.327];
    %     C = eye(25);
    % elseif model == "model4"
    %     %C = [0 0 -0.098 0.269;0 0 0.080 0.327];
    %     C = eye(20);
    % elseif model == "model5"
    %     %C = [0 0 -0.098 0.269;0 0 0.080 0.327];
    %     C = eye(40);
    % elseif model == "model6"
    %     %C = [0 0 -0.098 0.269;0 0 0.080 0.327];
    %     C = eye(60);
    % end
    a = csvread(input_folder + "/a.csv"); %A_matrix
    C = eye(size(a,1));
    b = csvread(input_folder + "/b.csv"); %B-vector
    x0 = csvread(input_folder + "/x0.csv"); %X0 values
end

