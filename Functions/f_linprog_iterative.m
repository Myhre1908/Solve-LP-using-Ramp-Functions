
function [xsave,usave,time] = f_linprog_iterative(c,H,h,G,Su,x0,A,B,tend)

    options = optimoptions('linprog', 'Display', 'off'); 
    nx = size(A,2); nu = size(B,2);
    dim_x = 2*nx;
    dim_u = 2*nu;
    xsave = zeros(nx,tend);
    usave = zeros(nu,tend);
    xsave(:,1) = x0;
    ll = 0;
    
    tic
    for pp = 2:tend+1
        try
            [z, fval, exitflag, output] = linprog(c, H, h, G, Su*x0, [], [], options);
            if isempty(z) || exitflag <= 0
                ll = 1;  % Set your condition flag
                %disp('linprog did not return a valid solution or failed');
                break;   % Exit the loop or take other actions
            end
        catch ME
            ll = 1;  
            disp('An error occurred during linprog execution:');
            disp(ME.message);  % Print the error message
            disp('Error identifier:');
            disp(ME.identifier);  % Print the error identifier
            disp('Stack trace:');
            disp(ME.stack);  % Display where the error occurred
            break;   
        end
        z_values = z(1:dim_x/2,:) - z(dim_x/2+1:dim_x,:);
        u_values = z(dim_x+1:dim_x+dim_u/2,:) - z(dim_x+1+dim_u/2:end,:);
        u = zeros(nu,1);
        x = zeros(nx,1);
        
        for kk = 1:nx
            x(kk,:) = z_values(kk:nx:end);
        end
        for ff = 1:nu
            u(ff,:) = u_values(ff:nu:end);
        end
        x0 = x;
        u = u_values(:,end);
        xsave(:,pp) = x0;
        usave(:,pp-1) = u;
    end
    if ll == 0
        time = toc;
    else
        time = NaN;
    end
end