function mpc = mpc_ingredients(A,B,Hx,hx,Hu,hu,CIS_H,CIS_h,x_ref_lin,u_ref_lin,Q,R,N)
    n = size(A,2);
    m = size(B,2);
    n_c_u = size(Hu,1);
    n_c_x = size(Hx,1);
    
    % Costruzioni delle matrici aumentate
    A_cal = zeros((N+1)*n,n);
    for i=0:N
        A_cal(i*n+1:(i+1)*n,:) = A^i;
    end
    
    B_cal = zeros((N+1)*n,N*m);
    for i=1:N
        for j=i-1:-1:0
            B_cal(i*n+1:(i+1)*n,(i-1-j)*m+1:(i-j)*m) = (A^j)*B;
        end
    end
    
    Hu_tilde = zeros(n_c_u*N,N*m);
    for i=0:N-1
        Hu_tilde(i*n_c_u+1:(i+1)*n_c_u,i*m+1:(i+1)*m) = Hu;
    end
    
    hu_tilde = zeros(n_c_u*N,1);
    for i=0:N-1
        hu_tilde(i*n_c_u+1:(i+1)*n_c_u,1) = hu;
    end
    
    Hx_tilde = zeros(n_c_x*(N+1),n*(N+1));
    for i=0:N
        Hx_tilde(i*n_c_x+1:(i+1)*n_c_x,i*n+1:(i+1)*n) = Hx;
    end
    
    C_H_tilde = zeros(size(CIS_H,1),n*(N+1));
    C_H_tilde(1:end,end-n+1:end) = CIS_H;
    
    Hx_tilde = [Hx_tilde; C_H_tilde];
    
    hx_tilde = zeros(n_c_x*(N+1),1);
    for i=0:N
        hx_tilde(i*n_c_x+1:(i+1)*n_c_x,1) = hx;
    end
    
    hx_tilde = [hx_tilde;CIS_h];
    
    [~,S,~] = dlqr(A,B,Q,R);
    
    Q_tilde = zeros((N+1)*n,(N+1)*n);
    for i=0:N
        Q_tilde(i*n +1:(i+1)*n,i*n +1:(i+1)*n) = Q;
    end
    Q_tilde(N*n +1:(N+1)*n,N*n +1:(N+1)*n) = S;
    
    R_tilde = zeros(N*m,N*m);
    for i=0:N-1
        R_tilde(i*m + 1:(i+1)*m,i*m +1:(i+1)*m) = R;
    end
    
    mpc = struct;
    
    % Matrici per la funzione di costo
    mpc.F = R_tilde + B_cal' * Q_tilde * B_cal;
    mpc.f_base = B_cal' * Q_tilde' * A_cal;
    
    % Matrici per i vincoli
    mpc.A_ineq = [Hu_tilde;Hx_tilde*B_cal];
    mpc.b_ineq_base = [hu_tilde;hx_tilde];
    mpc.b_ineq_x0_factor = [zeros(n_c_u*N,n);Hx_tilde*A_cal];
end