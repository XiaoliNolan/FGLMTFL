function [x] = FLADMM(A,y,s1,s2,D,opts)

    Print = opts.Print;
    Tasks = opts.Tasks;
    Dimen = opts.Dimen;
    MaxIter = opts.MaxIter;
    rho = opts.rho;
    Tol = opts.Tol;
    ind = opts.ind;
    
    Z = D*D;
    dZ = diag(Z);
    ZZ = -Z + diag(dZ);
    
    I = eye(Dimen); 
    bv = A'*y;
    for t = 1:Tasks
        CC = A'*A + rho*(1+dZ(t))*I;
        Ch_R{1,t} = dpotrf_(CC);
    end
    
    pk = zeros(Dimen,Tasks);
    v = zeros(Dimen,Tasks);
    u = zeros(Dimen,Tasks);
    q = zeros(Dimen,Tasks);
    x = zeros(Dimen,Tasks);
    x_p = zeros(Dimen,Tasks);
    
    ind_work(1:2,:)=ind(1:2,:);
    ind_work(3,:)=ind(3,:) * (s2 / rho);
    
    if Print
        fprintf('%3s\t%10s\t%10s\t%10s\t%10s\n','iter','funv', 'Tol tv', 'td', 'Tol tz');
    end
      
    for iter = 1:MaxIter
        % x-update
        TTv = v*D;
        TTpk = rho*pk*D;
        TTx = rho*x*ZZ; 
        xt = bv-u+rho*q+TTpk-TTv+TTx;
        for t = 1:Tasks
            x(:,t) = dpotrs_(Ch_R{1,t},xt(:,t));
        end
                    
        q = x + u/rho; 
        q = shrinkage_21(q,s1/rho);
        
        Tx = x*D;
        pk = Tx + v/rho;   
        pk = shrinkage_g(pk, ind_work, size(ind_work,2), Dimen, Tasks);
        
        u = u + rho*(x-q); 
        v = v + rho*(Tx-pk);   
        
        funv(iter) = norm(x-q);
        td = norm(x-x_p);
        
        if iter > 1
            tv = max(1,funv(iter-1));
            tz = max(1,td);
            if Print
                fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n',iter,funv(iter), Tol*tv,td,Tol*tz);
             end
            if funv(iter) < Tol*tv && td < Tol*tz
                break;
            end
        end
        
        x_p = x;
        
    end
    x = q;
             
function y = shrinkage(a, kappa)
    y = max(0, a-kappa) - max(0, -a-kappa);
    
function w = shrinkage_21(W,lambda_3)

    for i = 1:size(W,1)
        w_1 = W(i,:);
        nm = norm(w_1, 2);
        if nm == 0
            w_2 = zeros(size(w_1));
        else
            w_2 = max(nm - lambda_3, 0)/nm * w_1;
        end
        w(i,:) = w_2;
    end

function x = shrinkage_g(v, ind, nodes, nd, nt)
    x = reshape(v,nd*nt,1);
    for i = 1:nodes
        temp = ind(:,i);
        j = temp(1):temp(2);
        twoNorm = sqrt(sum(x(j,1).^2));
        lambda = temp(3);
        if twoNorm > lambda
            ratio = (twoNorm - lambda)/twoNorm;
            x(j,1) = x(j,1) * ratio;
        else
            x(j,1) = 0;
        end
    end
    x = reshape(x,nd,nt);


    