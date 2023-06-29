function [G,g] = cis(A,B,x_ref,u_ref,Fx,fx,Fu,fu,Q,R)

% Controllore LQR
K = -dlqr(A,B,Q,R); % Segno - per usare con + sotto

% Matrice A del sistema controllato con LQR
A_lqr = A+B*K;

% Vincoli su stato e ingresso
F = [Fx; Fu*K];
f = [fx; fu + Fu*(K*x_ref - u_ref)];

% Calcolo CIS (G*x<=g)
% Inizializzazione
CIS_poly_prev = Polyhedron(); % Poliedro vuoto
CIS_poly_curr = Polyhedron(F,f); % Poliedro fornito dai vincoli appena trovati
% La prima condizione del ciclo while (CIS_poly_prev.isEmptySet) serve per 
% non uscire alla prima iterazione
while CIS_poly_prev.isEmptySet || CIS_poly_prev ~= CIS_poly_curr
    % memorizzo il vecchio set candidato
    CIS_poly_prev = CIS_poly_curr;
    % Calcolo il nuovo set candidato (G_hat*x<=g_hat)
    G_hat = [CIS_poly_curr.A * A_lqr; F]; % tutti i vincoli stato
    g_hat = [CIS_poly_prev.b + CIS_poly_curr.A*B*(K*x_ref - u_ref); f];
    CIS_poly_curr = Polyhedron(G_hat,g_hat); % Poliedro con nuovi vincoli
end

% Disequazioni che descrivono il CIS
% Usiamo queste funzioni di MPT3 e non le implementiamo noi da zero perchè 
% queste tolgono in automatico i vincoli ridondanti, altrimenti il calcolo 
% diventa molto più oneroso
G = CIS_poly_curr.A;
g = CIS_poly_curr.b;

end

