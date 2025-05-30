clearvars
clc
addpath("Functions");

model = "model3";
% model = "model3";

Q_epsilon = 1e-1;
% Q_epsilon = 3; %Uncomment this to show the affect of the scaling
%algorith
c_scale = 10; %Myhres alpha value

[Hess,Hess_inv,G_eq,g_eq,h_ineq,H_ineq,C_A,A,B,C,x00,nx,nu,ny] = set_up_matrices(model,Q_epsilon,1);
x0 = x00(1,:)'; %Initial values from test bank
c = generate_random_linear_weight(nx,nu);

c = c*c_scale;
tend = 15;
max_iter = 400;

[xsave,usave,time2]               = f_linprog_iterative(c,H_ineq,h_ineq,G_eq,C_A,x0,A,B,tend);
[xsave2,usave2,time3,k3]          = f_run_quadprog_iterative(G_eq,C_A,H_ineq,h_ineq,c,Hess,A,B,x0,tend,max_iter);
[x,u,time1,k]                     = f_Run_Ramp_LP(G_eq,C_A,H_ineq,h_ineq,c,Hess_inv,Hess,A,B,x0,tend,max_iter,Q_epsilon);
compare_solutions(x,xsave,xsave2,time1,time2,time3)

disp(k) %Iteration for each z_k for the ramp function solver
disp(k3) %Iteration for each z_k for quadprog - Used to test Myhre's Algorithm 1 and 2
% 
% f_plot(x,nx,0," ")
% f_plot(xsave,nx,1,"--")

%%
function f_plot(x,nx,p,line)
    if p == 0
        hold off
    end
    for i = 1:nx
        if line ~= "--"
            plot(x(i,:))
        else
            plot(x(i,:),"--")
        end
        hold on
        title("State trajectories")
    end
end

function c = generate_random_linear_weight(nx,nu)
    c = zeros(2*nx*2*nu:1,1);
    c(1:nx) = (1+(10-1) * rand(nx,1));
    c(nx+1:2*nx) = c(1:nx);
    c(2*nx+1:2*nx+nu) =rand(nu,1);
    c(2*nx+nu+1:2*nx+2*nu) = c(2*nx+1:2*nx+nu);
end

function compare_solutions(x,xsave,xsave2,time1,time2,time3)
    difference = abs(x-xsave);
    difference3 = abs(x-xsave2);
    difference_mathworks = abs(xsave2-xsave);
    
    max_diff1 = max(difference(:));
    max_diff3 = max(difference3(:));
    max_diff4 = max(difference_mathworks(:));
    disp("Max differences:")
    disp("Ramp and linprog    :" + max_diff1)
    disp("Ramp and quadprog   :" + max_diff3)
    disp("quadprog and linprog:" + max_diff4)
    disp(" ")
    disp("Ramp solver time [s]: " + time1)
    disp("linprog time     [s]: " + time2)
    disp("quadprog time    [s]: " + time3)
    disp(" ")
end

