% =========================================================================
%  STOCHASTIC SIMPLICIAL SIR CONTAGION MODEL (S3IR-CM)
%  SIR Dynamics on a Simplicial Complex with Higher-Order Interactions
%  Using emperical data sets
 % =========================================================================
%
%  Model Description:
%  ------------------
%  This script simulates contagion dynamics, incorporating pairwise (D=1),
%  triangular (D=2), and tetrahedral (D=3) interactions. Each dimension
%  adds a higher-order contagion mechanism beyond standard pairwise contact.
%
%  State Variables (per node i, per time step t):
%  -----------------------------------------------
%    S(i,t) = 1  if node i is Susceptible   I(i,t)=0, R(i,t)=0
%    I(i,t) = 1  if node i is Infected     
%    R(i,t) = 1  if node i is Recovered   
%
%  Infection Mechanism (Poisson process):
%  --------------------------------------
%    λᵢ infection probability updates at each time step if node is part of
%    any interaction
%
%    p_infect = 1 - exp(-λᵢ*dt) 
%
%  Parameters:
%  -----------
%    β1, β2, β3  : infection rates for 1-, 2-, 3-simplex respectively (per unit time)
%    γ            : recovery rate (per unit time)
%    dt           : time step 
%    σ1,σ2,σ3: noise parameters

clc
rng(100,"twister")

%% ===========================================================
% SECTION 1: Import data
%=============================================================
% The simplicial complex is defined by three files containg edges, triangles and tetrahedra respectively:
%    edges_invs.csv       each row is a pair   [i, j]
%    triangles_invs.csv    each row is a triple [i, j, k]
%    Tetrahedra_invs.csv   each row is a quad   [i, j, k, l]
%
%    n = total number of nodes, inferred from the maximum node index
%=============================================================
edges=readmatrix('edges_invs.csv');              
triangles = readmatrix('triangles_invs.csv');   
tetrahedron= readmatrix('tetrahedra_invs.csv');
n = max([edges(:); triangles(:); tetrahedron(:)]);

%% ========================================================================
%  SECTION 2: Simulation time settings
%  ========================================================================
%  T  : total simulation time (continuous-time units)
%  dt : discrete time step
%  N  : number of discrete time steps = T/dt
% =========================================================================
T=10;                   % total simulation time
dt=0.001;               % stepsize
N=T/dt;                 % total discrete time steps

%% =======================================================================
%SECTION 3: Epidemic parameters
%  ========================================================================

bEpidemic = false;
if bEpidemic    
    fname = 's3ir_p';
else            
    fname = 's3ir_e'; 
end
if bEpidemic   % outbreak will spread
 gamma=0.5;    % recovery rate
 beta_1=0.1;   % Pairwise infection rate
 beta_2=0.2 ;  %tri angular infection rate
 beta_3=0.3;   %tetrahedral infection rate 
 sigma1=0.03; sigma2=0.02;  sigma3=0.01; % noise parameters
else           % outbreak will die out
 gamma=0.5; beta_1=0.01;beta_2=0.02 ;beta_3=0.03;  
 sigma1=0.01; sigma2=0.02;  sigma3=0.01;
end
%% ========================================================================
%  SECTION 4: Precompute neighborhood list
%  ========================================================================
%  For efficiency, we pre-compute for each node i the list of simplex rows
%  (edges, triangles, tetrahedra) that contain it. This avoids calling
%  find(any(...)) inside every time step.
%
%  edge_find{i}  = row indices in 'edges'       that contain node i
%  tri_find{i}   = row indices in 'triangles'   that contain node i
%  tetra_find{i} = row indices in 'tetrahedron' that contain node i
% =========================================================================
edge_find  = cell(n, 1);
for i = 1:n
    edge_find{i}  = find(any(edges == i, 2));
end

tri_find   = cell(n, 1);
for i = 1:n
    tri_find{i}   = find(any(triangles == i, 2));
end

tetra_find = cell(n, 1);
for i = 1:n
    tetra_find{i} = find(any(tetrahedron == i, 2));
end
seed_nodes = randperm(n, 2);                      % initial infection
%% ========================================================================
%  SECTION 5 : SIMULATION  — 3-simplex
%  ========================================================================
%  Run the loop for multiple realizations.
%  All connections in a 3-simplex contribute to infection.
%  pairwise interactions contribute to infection if 
%  a susceptible node i is exposed if any of its edge-neighbours is infected.
%  A triangle is "fully active" for node i if both other triangle members
%  are simultaneously infected
%  A tetrahedron is "fully active" for node i if all three other members
%  are simultaneously infected 
%

NR=25;                                           
f = figure('visible', 'on');
hold on
for run=1:NR                                      % multiple iterations

I = zeros(n,N);                                   % infection state matrix
R = zeros(n,N);                                   % recovery state matrix

I(seed_nodes, 1) = 1;                             % Set all seed node to 1 initially (infected)
R(:,1) = 0;                                       % no node is recovered at t=0

for t = 1:N-1
    I_new = I(:, t);                              % initiallize next state using current state
    R_new = R(:, t);
    xi_t  = randn;                                  
    beta1 =  beta_1 + sigma1*xi_t;          % perterbued infection rates
    beta2 =  beta_2 + sigma2*xi_t;
    beta3 =  beta_3 + sigma3*xi_t;

    for i = 1:n
        if I(i,t) == 0                              % susceptible node
            lambda_i = 0;
            % STEP 1: Edge-based (pairwise) contagion transmission
            %   For each edge containing node i, check if the neighbour
            %   is infected. Each infected edge-neighbour contributes β₁
            %   to the infection rate λᵢ.

            for e = edge_find{i}'
                pair = edges(e, :);
                nhbr = pair(pair ~= i);
                if I(nhbr, t) == 1 
                    lambda_i = lambda_i + beta1;
                 end
            end
            %   STEP2: Triangle-based infection
            %   For triangular face containing node i, check if the two
            %   neighbours are infected, the triangle is fully active for node i and contribute at rate beta2.

            lambda_tri = 0;
            for k = tri_find{i}'
                tri    = triangles(k, :); 
                others = tri(tri ~= i);
                if sum(I(others, t)) == 2 
                    lambda_tri = lambda_tri + 1;
                end
            end
            lambda_i = lambda_i + beta2 * lambda_tri;       % update infection probability

            % STEP3: Tetrahedron-based infection
            % The tetrahedron is fully active for node i if three
            % neighbours are infected. this is the maximum contagion
            % pressure.
            lambda_tet = 0;
            for h = tetra_find{i}'
                tet    = tetrahedron(h, :);
                others = tet(tet ~= i);
                if sum(I(others, t)) >= 2 
                    lambda_tet = lambda_tet + 1;
                end
            end
                                                         
            lambda_i = lambda_i + beta3 * lambda_tet;     % update infection probability
            p_infect = 1 - exp(-lambda_i*dt);
            if rand < p_infect
                I_new(i) = 1;
                R_new(i) = 0;
            end

        elseif I(i,t) == 1                               % If the node is infected
            % STEP4: Recovery
            if rand < gamma*dt
                I_new(i) = 0;
                R_new(i) = 1;
            end
        end
    end
    % STEP5: Update states
    I(:, t+1) = I_new;
    R(:, t+1) = R_new;
end

I_t = mean(I, 1);                                       % Find avarage infection
R_t = mean(R, 1); 
S_t = 1 - I_t - R_t;

%%=========================================================================
% SECTION 6: plot the graph
%===========================================================================

plot(dt:dt:T,I_t);
ylim([0 1])
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1]; 
xlabel('$t$','Interpreter','latex', 'FontSize',15);
title('$I(t)$', 'Interpreter','latex', 'FontSize',15);
grid on;
box on;
exportgraphics(f, 'inf_ext_invs.png');
end
hold off