function x_dot = cstr(t,x,u,q0,V,k0,ER,dHr,UA, ro, Cp, tau_c, CAf, Tf)
%CSTR Modello nello spazio degli stati di un reattore CSTR

% Equazioni di stato
x_dot = zeros(3,1);
x_dot(1) = (q0/V)*(CAf - x(1))-k0*exp(-ER/x(2))*x(1);
x_dot(2) = (q0/V)*(Tf - x(2)) + (((-dHr)*k0)/(ro*Cp))*exp(-ER/x(2))*x(1) + UA/(V*ro*Cp)*(x(3) - x(2));
x_dot(3) = (u - x(3))/tau_c;
end