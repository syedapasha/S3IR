%% ═══════════════════════════════════════════════════════════════════════
%  Stochastic SIR Model on Simplicial Complex (Real Dataset)
%__________________________________________________________________________
%Description:
%  This script simulates contagion dynamics, incorporating pairwise (D=1),
%  triangular (D=2), and tetrahedral (D=3) interactions in a real dataset
%  using Euler Maryama scheme.(Infection rates are updated at D=1,2,3 respectively 
%  and then total force of infection is used to simulate Infection dynamics 
% corresponding to each node using EM scheme)
%__________________________________________________________________________
%   Higher-order interactions:
%    - Pairwise  (edges)        → beta1
%    - Triangle  (2-simplex)    → beta2
%    - Tetrahedron (3-simplex)  → beta3
%
%  Input CSV files:
%   The simplicial complex is defined by three files containg edges, triangles and tetrahedra respectively:
%    edges_invs.csv       each row is a pair   [i, j]
%    triangles_invs.csv    each row is a triple [i, j, k]
%    Tetrahedra_invs.csv   each row is a quad   [i, j, k, l]
%    n = total number of nodes, inferred from the maximum node index
%  
%% ═══════════════════════════════════════════════════════════════════════

clc; clear; close all;
rng(100,"twister"); 

%% ── 1. LOAD CSV DATA ────────────────────────────────────────────────────


edges      = readmatrix('edges_invs.csv');        
triangles  = readmatrix('triangles_invs.csv');    
tetrahedra = readmatrix('tetrahedra_invs.csv');   

% Total number of nodes
n = max([edges(:); triangles(:); tetrahedra(:)]);

%% ── 2. BUILD LOOKUP STRUCTURES ──────────────────────────────────────────

neighbors  = cell(n, 1);   % neighbors{i}  = list of pairwise neighbors
tri_find   = cell(n, 1);   % tri_find{i}   = list of triangle indices containing i
tetra_find = cell(n, 1);   % tetra_find{i} = list of tetrahedron indices containing i

% Pairwise neighbors from edge list
for e = 1:size(edges, 1)
    u = edges(e,1);  v = edges(e,2);
    neighbors{u}(end+1) = v;
    neighbors{v}(end+1) = u;
end

% Triangle membership
for t = 1:size(triangles, 1)
    for k = 1:3
        tri_find{triangles(t,k)}(end+1) = t;
    end
end

% Tetrahedron membership
for h = 1:size(tetrahedra, 1)
    for k = 1:4
        tetra_find{tetrahedra(h,k)}(end+1) = h;
    end
end

%% ── 3. AVERAGE NODE TO FACET DEGREE──────────────────────────────────────────
% pre computed values
k_avg     = 16.413;
k_tri_avg = 1.7283;
k_tet_avg = 1.5217;




%% ── 4. MODEL PARAMETERS ─────────────────────────────────────────────────
bEpidemic=true;
if bEpidemic
beta1  = 0.05;    % pairwise transmission rate
beta2  = 0.1;    % triangle transmission rate
beta3  = 0.15;    % tetrahedron transmission rate
gamma     = 0.001;    % recovery rate  
sigma1  = 0.01; sigma2=0.01; sigma3=0.01;
else
beta1  = 0.005;    % pairwise transmission rate
beta2  = 0.01;    % triangle transmission rate
beta3  = 0.015;    % tetrahedron transmission rate
gamma     = 0.01;    % recovery rate  
sigma1  = 0.1; sigma2=0.1; sigma3=0.1;    % stochastic noise intensity
end
dt     = 0.001;    % Euler-Maruyama time step
T  = 100;      % total simulation time
N = round(T / dt); %number of time steps 
time   = (0:N-1) * dt;



%% ── 5. INITIAL CONDITIONS ───────────────────────────────────────────────
                    
I   = zeros(n, N);   % I(i,t) = infection probability of node i at step t
R   =zeros(n, N); 
S   =zeros(n, N);
% I0_frac = 0.03;              % seed 3 nodes as infected
% 
% seed_nodes = randperm(n, max(1, round(I0_frac * n)));
% I(seed_nodes, 1) = 1;


%% ── 6. EULER-MARUYAMA SDE SIMULATION ───────────────────────────────────

I(:,1)=0.3;
R(:,1)=0;
S(:,1)=1-I(:,1)-R(:,1);
for t = 1:N-1

    I_t = I(:, t);   % state at current step(current states)
    R_t = R(:, t);
    S_t = S(:,t);
     I_new = I_t;       % preallocate next step(new states)
     R_new=R_t;
     S_new=S_t;
    for i = 1:n

        S_t(i) = 1 - I_t(i)-R_t(i);   % susceptible fraction at node i

        % ── Pairwise force of infection ──────────────────────────────
        if ~isempty(neighbors{i})
            lambda_pair = (beta1 / k_avg) * sum(I_t(neighbors{i}));
        else
            lambda_pair = 0;
        end

        % ── Triangle infection ─────────────────────────────────────────────
        n_active_tri = 0;
        for h = tri_find{i}
            others = triangles(h, triangles(h,:) ~= i);
            if sum(I_t(others)) == 2       %  2 others infected
                n_active_tri = n_active_tri + 1;
            end
        end
        lambda_tri = (beta2 /k_tri_avg) * n_active_tri;

        % ── Tetrahedron infection ──────────────────────────────────────────
        n_active_tet = 0;
        for h = tetra_find{i}
            others = tetrahedra(h, tetrahedra(h,:) ~= i);
            if sum(I_t(others)) == 3       % majority infected
                n_active_tet = n_active_tet + 1;
            end
        end
        % infection update
        lambda_tet = (beta3 / k_tet_avg) * n_active_tet;

        % ── Total force of infection ──────────────────────────────────
       lambda_i = lambda_pair + lambda_tri + lambda_tet;
    

        % ── Euler-Maruyama update ─────────────────────────────────────
        

        dW     = sqrt(dt) * randn;
        
        diff=S_t(i)*(sigma1*I_t(i)+sigma2*I_t(i)^2+sigma3*I_t(i)^3)*dW;
        drift  = (lambda_i * S_t(i)  -  gamma * I_t(i)) * dt;
       
        I_new(i) =  I_t(i) + drift + diff;
     
        R_new(i) = R_t(i)+ gamma * I_t(i) * dt;
        S_new(i)=1-I_new(i)-R_new(i);
    end

    I(:, t+1) = I_new;
    R(:, t+1) = R_new;
    S(:, t+1) = S_new;
   
end
%% ── 7. AGGREGATE RESULTS ────────────────────────────────────────────────

I_mean = mean(I, 1);     % mean prevalence across all nodes
R_mean= mean (R,1);
S_mean= mean (S,1);


%% 8. PLOTS ────────────────────────────────────────────────────────────

% ── Figure 1: susceptible nodes ─────────────────────────────
figure(1); clf;
plot(time,S); axis([0 T 0 1]);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1]; 
xlabel('$t$','Interpreter','latex', 'FontSize',15);
title('$S(t)$', 'Interpreter','latex', 'FontSize',15);
grid on; box on;
%exportgraphics(figure(1), 's3irCM_S_invs.png');
% ── Figure 2: Per-node infection ─────────────────────────────
figure(2); clf;
plot(time,I); axis([0 T 0 1]);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1]; 
%plot((1:N)*dt,R); axis([0 T 0 1]);
xlabel('$t$','Interpreter','latex', 'FontSize',15);
title('$I(t)$', 'Interpreter','latex', 'FontSize',15);
grid on; box on;
%exportgraphics(figure(2), 's3irCM_I_invs.png');
% ── Figure 3: recovered nodes─────────────────────────────
figure(3); clf;
plot(time,R); axis([0 T 0 1]);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1]; 
%plot((1:N)*dt,R); axis([0 T 0 1]);
xlabel('$t$','Interpreter','latex', 'FontSize',15);
title('$R(t)$', 'Interpreter','latex', 'FontSize',15);
grid on; box on;
%exportgraphics(figure(3), 's3irCM_R_invs.png');
% ── Figure 4: Mean prevalence ───────────────────────────
figure(4); clf;
hold on;
plot(time, I_mean, 'r-',  'LineWidth', 1.0, 'DisplayName', 'Mean I(t)');
plot(time, R_mean,  'g-', 'LineWidth', 1.0, 'DisplayName', 'Mean R(t)');
plot(time, S_mean,  'b-', 'LineWidth', 1.0, 'DisplayName', 'Mean S(t)');
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1]; 
xlabel('$t$','Interpreter','latex', 'FontSize',15);
title('Mean prevalence', 'Interpreter','latex', 'FontSize',15);
legend('show', 'Location', 'east');
grid on; box on;
exportgraphics(figure(4), 's3irCM_M_invs.png');