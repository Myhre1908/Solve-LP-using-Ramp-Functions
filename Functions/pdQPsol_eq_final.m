function [v,lam,mu,Qmat0i,y]=pdQPsol_eq_final(IGis,g,h,GQGi,HQG,QG,QH,nc,ic,ik,q_epsilon)
%Ramp-based QP solution for (strictly) positive definite Hessians including
%equality constraints
% J = 0.5z'*Hess*z
%s.t. Gz = g
%     Hz <=h
%
%Input:
%IGis = eye(nc)+H*Hess^{-1}*H'-H*Hess^{-1}*G'*(G*Hess^{-1}*G')^{-1}*H*Hess^{-1}*H';
%nc = number of constraints (row dimension of IGis)
%...
%GQGi = inv(G*Hess^{-1}*G')
%HQG = H*Hess^{-1}*G'
%QG = Hess^{-1}*G'
%QH = Hess^{-1}*H'

%Output
%z   Optimal value of z
%lam    Lagrange multipliers for inequaliies
%mu     Lagrange multipliers for equalities
%flag = 0 Optimal solution found
%flag = -1  Infeasible
%flag = -2  Appears to have entered indefinite loop


rank2 = false;
solved = false;
flag = 0;

if nc < 100
    prioritize_add_elements = true;
else
    prioritize_add_elements = false;
end

ix = 0;
qC = [-1 0;0 1];
qV0 = zeros(2,nc);

if ic == 1
    y0i = HQG*GQGi*g-h; %m4 for g ~= 0
else
    y0i = -h; %if g (g-Gzc) = 0 then m5 becomes -h
end

min_treshold = 1e5*eps; 
y0 = y0i;
actset = zeros(nc,1);
val = 1;
Qmat0i = speye(nc); 
Min_qdiv = 1e-4;
% Min_qdiv = q_epsilon;

while (~solved)
    ix = ix+1;
    if (ix > 8*nc) %Add a pertubation to the algorithm
        if nc < 100 && qc == -1 %Want it to add active sets twice
            prioritize_add_elements = false;
        else
            prioritize_add_elements = true;
        end
    end

    if (ix > (8+val)*nc)%Reverts the pertubation so Q doesnt become singular
        if nc < 100
            prioritize_add_elements = true;
        else
            prioritize_add_elements = false;
        end
    end
   
    if (ix > 11*nc)
       flag = -2;
       break
    end
    if (ix == 1)
        y = y0;
    elseif ~rank2
        y = y0-y0(iz)*vAd;
        y0 = y;
    else
       y = Qmat0i*y0i;
       y0=y;
    end
    lam = y.*actset;

    if prioritize_add_elements %Try to not be confused by this :)
        %Maybe start from ~prioritize_add_active_elements, as this is the
        %old algorithm
        if ix < 8*nc
            [i2,iz] = max(y-lam);
        else
            chance = rand(1);
            [y_sorted, idx] = sort(y - lam, 'descend'); % Sort in descending order
            if chance > 2/3
                i2 = y_sorted(2);  % Get the second highest value
                iz = idx(2);       % Get its corresponding index
            else
                i2 = y_sorted(1);  % Get the highest value
                iz = idx(1);       % Get its corresponding index
            end
        end

        if i2 < min_treshold
            i2 = [];
        end

        if (~isempty(i2)) 
            actset0 = actset;
            actset(iz) = 1;
            qc = -1;
        elseif min(lam) < -min_treshold %&& ~rank2_activated
            [~,iz] = min(lam);
            actset(iz) = 0;
            qc = 1; % Remove constraint from active set
        else  
             iz = [];
        end



      elseif ~prioritize_add_elements
        [i1,i1z] = min(lam);
        if i1 < -min_treshold  
          iz  = i1z;
          actset0 = actset;
          actset(iz) = 0;
          qc = 1; 

        else
            [i2,i2z] = max(y-lam);
            if (i2 <= min_treshold)
                i2 = [];
            end
            if(~isempty(i2))
                iz  = i2z;
                actset0 = actset;
                actset(iz) = 1;
                qc = -1; 
            else
                iz = [];
            end
        end
    end


    if(~isempty(iz))
        qu = IGis(:,iz);
        vA = Qmat0i*qu;
        qdiv = qc+vA(iz);
        if abs(qdiv) < Min_qdiv
            rank2 = true;
        else
            rank2 = false;
        end

        if ~rank2 
            vAd = (1/qdiv)*vA;
            Qmat0i = sparse(Qmat0i-vAd * Qmat0i(iz, :));
        else %rank2
            actix = find(actset0==1);
            dy = Qmat0i*qu;
            dy = clean(dy,1e-6);
            dyn = find(dy(actix) < 0);
            if (isempty(dyn))
                flag = -1;
                break
            end
            values = abs(y(actix(dyn)) ./ dy(actix(dyn)));
            [~, sortedIndices] = sort(values, 'ascend');        
            if actix(dyn(sortedIndices(1))) == iz
                iz3 = sortedIndices(2);
            else
                iz3 = sortedIndices(1);
            end
            iz4 = actix(dyn(iz3));

            if iz4 == iz %If iz4 = iz => reassign iz4 to the second most extreme value
                iz3 = sortedIndices(2);
                iz4 = actix(dyn(iz3));
            end
            actset(iz4)=0;
            qu2 = IGis(:,iz4);
            qU = [qu qu2];
            qV = qV0;
            qV(1,iz) = 1;
            qV(2,iz4) = 1;
            % Qmat0i = Qmat0i * qU * 1/(qC + Qmat0i([iz, iz4],:)*qU) * Qmat0i([iz, iz4],:);
            Qmat0i = sparse(Qmat0i - Qmat0i * qU * 1/(qC + qV*Qmat0i*qU) * qV*Qmat0i);  
        end
    else
        solved = true;
    end
        
 end
    if (flag == 0)
        mu = -GQGi*(HQG'*lam+g);
        v = -QG*mu-QH*lam;
    else
        if flag == -1
            disp("Infeasible")
        else
            disp("Infinity loop")
        end
        lam = [];
        mu = [];
        v = [];
    end
end

function x = clean(x, threshold)
    x(abs(x) < threshold) = 0;
end