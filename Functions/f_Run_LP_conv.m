function [xsave,usave,time,isave,z_save] = f_Run_LP_conv(G,Su,H,h,c,Hi,Hess,A,B,x0,tend,max_iter,Q_epsilon)
    n = 1;
    ny = size(A,2);
    nu = size(B,2);
    nc = size(H,1);
    nv = n*(ny*2+nu*2);
    dim_x = 2*ny*n;
    dim_u = 2*nu*n;
    xsave = zeros(ny,tend+1);
    usave = zeros(nu,tend);
    isave = zeros(1,tend);
    xsave(:,1)=x0;
    z_save = zeros(nv,tend);
    

    
    [GQGi,HQG,HQ,QH,QG,IGi] = create_matrices(Hi,G,H,nc);
    %IGi = M4
    M = QG*GQGi*HQG'-QH; %Myhre's M_\gamma

    h_HQc = h + HQ*c; 
    Hi_c = Hi*c;
    gv1 = zeros(size(x0));

    tic
    for ik = 1:tend
        zc = zeros(nv,1);
        % zcc = zeros(nv,1);
        ic = 0;  
        g = Su*x0;
        scale_zc = false;
        while true
            ic = ic + 1;
            if ic == 1
                %For first Iteration g = g - Gzc
                %However zc = 0 making it just g
                gv = g; 
                gv = clean(gv, 1e-6);
            else
                 %For remainding iterations g = g - Gzc
                 %However g-G*zc = 0 something that can be verified by
                 %testing
                gv = gv1;
            end
            if scale_zc %Myhre's Agortihm 2
                dy = y2-y;
                scale = find_optimal_lam_scaling(y, dy); 
                if scale < 5e3
                    zc = z2 + (z2-z)*scale;
                else
                    disp("Guessing we have found the correct active set")
                    break
                end
                hv = h_HQc - H*zc;
            else
                hv = h_HQc - H*zc;
            end 

            [v, lam, ~,Qmat0i,y] = pdQPsol_eq_final(IGi, gv, hv, GQGi, HQG, QG, QH,nc,ic,ik,Q_epsilon);

            if isempty(lam) %Ramp functions failed
                ik
                ic
                max(xsave(:,ik))
                break
            end
            z = v-Hi_c+zc;
            if max(Hess*v-c) < 1e-6
                break
            end

            if ic > 2 %Myhres Algorithm 1
                [z2,~,~,y2] = get_z(Qmat0i,h_HQc,H,z,M,Hi,c);
                diff = max(z-zc);
                diff2 = max(z2-z);
                y_signs_changes = check_signs_equal(y, y2);
                if y_signs_changes && abs(diff-diff2) < 1e-3
                    scale_zc = true;
                else
                    scale_zc = false;
                end
            end

            if ic == max_iter %So it doesnt enter an infinite loop
                break;
            end
            zc = z;
        end

        if ~isempty(lam)
            u =  z(dim_x+1:dim_x+dim_u/2,:) - z(dim_x+1+dim_u/2:end,:); %up - um
            x0 = z(1:dim_x/2,:) - z(dim_x/2+1:dim_x,:);

            xsave(:,ik + 1) = x0;
            usave(:,ik) = u;
            isave(1,ik) = ic;          
        else
            xsave = zeros(ny,tend+1);
            usave = zeros(nu,tend);
            isave = zeros(1,tend);
            break
        end
    end

    if ik == tend
        time = toc;
    else
        time = NaN;
    end

end

function x = clean(x, threshold)
    x(abs(x) < threshold) = 0;
end

function [GQGi,HQG,HQ,QH,QG,IGi] = create_matrices(Hi,G,H,nc)
    GQG = G*Hi*G'; GQGi = inv(GQG);
    HQH = H*Hi*H'; HQG = H*Hi*G';
    HQ = H*Hi; QH = HQ';
    GQ = G*Hi; QG = GQ';
    IGi = (eye(nc)-HQH+HQG*GQGi*HQG');%M4
end

function n = find_optimal_lam_scaling(y, dy)
    % Ensure dy is nonzero to avoid division by zero
    valid_indices = dy ~= 0;
    
    % Compute the minimum positive integer n that makes any element switch sign
    n_values = -y(valid_indices) ./ dy(valid_indices);
    
    % Select the smallest positive integer solution
    n = min(n_values(n_values > 0));
    % n = ceil(min(n_values(n_values > 0))); For integer
end

function [z,v,lam,y] = get_z(Qmat0i,h_HQc,H,z,M,Hi,c)
    lam = -Qmat0i*(h_HQc-H*z);
    y = lam;
    lam(lam<0) = 0;
    v = M*lam;
    z = v-Hi*c+z;
end

function isSameSign = check_signs_equal(y, y_possible_new)
    % Ensure vectors have the same length
    if length(y) ~= length(y_possible_new)
        error('Vectors must have the same length.');
    end
    
    % Check if all elements have the same sign
    isSameSign = all(sign(y) == sign(y_possible_new));
end
