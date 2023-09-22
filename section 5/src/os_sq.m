function [A,B,C,Q] = os_sq(U,kx,kz,R)
    n = length(U);
    k = sqrt(kx^2+kz^2);
    y = cheb_grid(n);
    D = cheb_diff(n);
    D2 = D^2;
    S = diag([0; 1 ./(1-y(2:n-1).^2); 0]);
    D4 = (diag(1-y.^2)*D^4 -8*diag(y)*D^3-12*D^2)*S;
    w = cheb_weights(n);

    U_tilde = U(2:n-1);
    D_tilde = D(2:n-1,2:n-1);
    D2_tilde = D2(2:n-1,2:n-1);
    D4_tilde = D4(2:n-1,2:n-1);
    w_tilde = w(2:n-1);
    Z = zeros(n-2,n-2);
    I = eye(n-2);
    Q = [(k^2*diag(w_tilde)+D_tilde'*diag(w_tilde)*D_tilde) Z;
          Z diag(w_tilde)]/(2*k^2);

    M = k^2*I-D2_tilde;
    M2 = k^4*I-2*k^2*D2_tilde+D4_tilde;
    invM = M\I;

    Aos = invM*(-1i*kx*diag(U_tilde)*M-1i*kx*diag(D2_tilde*U_tilde)-(1/R)*M2);
    Asq = -1i*kx*diag(U_tilde)-(1/R)*M;
    Ac = -1i*kz*diag(D_tilde*U_tilde);

    A = [Aos Z;
         Ac Asq];

    B = [1i*kx*invM*D_tilde invM*k^2 1i*kz*invM*D_tilde;
        1i*kz*I Z -1i*kx*I];

    C = [(1i*kx/k^2)*D_tilde (-1i*kz/k^2)*I;
        I Z;
        (1i*kz/k^2)*D_tilde (1i*kx/k^2)*I];
end