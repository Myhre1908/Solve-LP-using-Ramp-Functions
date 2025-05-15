function [x,usave,time,isave] = f_run_quadprog_iterative(G2,Su2,H2,h2,c2,Hess2,A2,B2,x0,tend,max_iter)
    n = 1;
    nx = size(A2,2); nu = size(B2,2);
    x00 = x0;
    dim = n*(2*nx+2*nu);
    options = optimoptions('quadprog', 'Display', 'off');
    dim_x = 2*nx*n;
    dim_u = 2*nu*n;
    Hess_inv = eye(size(Hess2))*1/Hess2(1,1);
    isave = zeros(1,tend);

    
    
    Inv_H = inv(Hess2);
    xsave = zeros(nx,tend);
    usave = zeros(nu,tend);
    
    tic
    for k = 1:tend
        z         = zeros(dim,max_iter);
        zc        = zeros(dim,1);
        gamma_old = zeros(dim,1);
        i = 1;
        ll = 0;
        
        while true
            var_change = (c2- Hess2'*zc);
            add_G = G2*Inv_H*var_change;
            add_H = H2*Inv_H*var_change;
            try
                [gamma, ~, exitflag, ~,lambdak] = quadprog(Hess2, [], H2, h2 + add_H, G2, Su2 * x0 + add_G, [], [], [], options);
                % if i == 3
                %     tols = 1e-6; % Set a tolerance level
                %     ineq_significant = lambdak.ineqlin(abs(lambdak.ineqlin) > tols);
                %     eq_significant = lambdak.eqlin(abs(lambdak.eqlin) > tols);
                %     disp(ineq_significant); % Only display multipliers above the threshold
                %     % disp(eq_significant);
                % end
                if isempty(gamma) || exitflag <= 0
                    ll = 1;  
                    break;   
                end
            catch
                ll = 1;  % Set your condition flag
                %disp('Error occurred during quadprog execution');
                break;   % Exit the loop or take other actions
            end
    
            z = gamma - Inv_H*var_change;
            zc = z;
            if (max(Hess2*gamma-c2) < 1e-8)
                % fprintf('Solved at iteration %d \n',i)
                break
            end
            
            if i == max_iter
                break
            end
            % max(gamma-gamma_old)
            gamma_old = gamma;
            i = i+1;
        end
        if ll == 0
            z_values = z(1:dim_x/2,:) - z(dim_x/2+1:dim_x,:);
            u_values = z(dim_x+1:dim_x+dim_u/2,:) - z(dim_x+1+dim_u/2:end,:);
            x0 = z_values;
        
            xsave(:,k) = x0;
            usave(:,k) = u_values;
            isave(k) = i;
        end
    end
    
    if ll == 0
        time = toc;
        % isave
    else
        time = NaN;
    end
    x = [x00,xsave];
end
