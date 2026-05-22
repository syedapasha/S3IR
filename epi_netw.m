% =========================================================================
%  SIR transmission on simplicial complex structure revelaed from  
%  emperical datasets
% =========================================================================
%
%  Model Description:
%  ------------------
%  This script simulates contagion dynamics, incorporating pairwise (D=1),
%  triangular (D=2), and tetrahedral (D=3) interactions, and higher-degree
%  interactions. Each dimension > 1
%  adds a higher-order contagion mechanism beyond standard pairwise contact
%
%  State Variables (node i, time step t):
%  -----------------------------------------------
%    S(i,t) = 1  if node i is Susceptible   I(i,t) = 0, R(i,t) = 0
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
%    βᵢ  : infection rates (for i-simplex) 
%    γ   : recovery rate
%    dt  : time step 
%    σᵢ  : noise parameters

clear all; close all; clc;
rng(100,"twister")

%% ===========================================================
% Import data
%=============================================================
% The simplicial complex is defined by three files containg edges, 
% triangles and tetrahedra respectively:
%
%    edges_invs.csv        each row is a pair   [i, j]
%    triangles_invs.csv    each row is a triple [i, j, k]
%    Tetrahedra_invs.csv   each row is a quad   [i, j, k, l]
%
%    n = number of nodes (inferred from the maximum node index)
%=============================================================
bSFHH = true;

if bSFHH    dname = 'sfhh';    else    dname = 'invs';  end

edges       = readmatrix(strcat('edges_',dname,'.csv'));              
triangles   = readmatrix(strcat('triangles_',dname,'.csv'));   
tetrahedron = readmatrix(strcat('tetrahedra_',dname,'.csv'));
five_face   = readmatrix(strcat('penta_',dname,'.csv'));
six_face    = readmatrix(strcat('hexa_',dname,'.csv'));
seven_face  = readmatrix(strcat('hepta_',dname,'.csv'));
    
n = max([edges(:); triangles(:); tetrahedron(:); five_face(:); six_face(:); seven_face(:)]);

%% ========================================================================
%  Simulation settings
%  ========================================================================
%  T  : simulation time 
%  dt : step size
%  N  : number of discrete time steps = T/dt
% =========================================================================
T = 100; 
dt= 1e-1; 
N = T/dt;  

%% ========================================================================
%  Epidemic parameters
%  ========================================================================

bEpidemic = false;

if bEpidemic    fname = 'Epi';   else    fname = 'noEpi';   end

bet_1 = 0.01;  bet_2 = 0.02;  bet_3 = 0.1; bet_4 = 0.1;  bet_5 = 0.1;  bet_6 = 0.1;
sgma1 = 0.01;  sgma2 = 0.01;  sgma3 = 0.01; sgma4 = 0.01;  sgma5 = 0.01;  sgma6 = 0.01; 

if bEpidemic            % outbreak will spread   
    gamma = 0.1;    
else                    % outbreak will die out
    gamma = 0.6; 
end

%% ========================================================================
%  Precompute neighborhood list
%  ========================================================================
%  For efficiency, we precompute for each node i the list of simplex rows
%  (edges, triangles, tetrahedra) that contain it. This avoids calling
%  find(any(...)) inside every time step.
%
%  edge_find{i}  = row indices in 'edges'       that contain node i
%  tri_find{i}   = row indices in 'triangles'   that contain node i
%  tetra_find{i} = row indices in 'tetrahedron' that contain node i
% =========================================================================

edge_find = cell(n,1);
for i = 1:n
    edge_find{i} = find(any(edges == i,2));
end

tri_find = cell(n,1);
for i = 1:n
    tri_find{i} = find(any(triangles == i,2));
end

tetra_find = cell(n,1);
for i = 1:n
    tetra_find{i} = find(any(tetrahedron == i,2));
end

penta_find = cell(n,1);
for i = 1:n
    penta_find{i} = find(any(five_face == i,2));
end

hexa_find = cell(n,1);
for i = 1:n
    hexa_find{i} = find(any(six_face == i,2));
end

hepta_find = cell(n,1);
for i = 1:n
    hepta_find{i} = find(any(seven_face == i,2));
end


% initial infection
if bSFHH    seed_nodes = randperm(n,10);    else    seed_nodes = randperm(n,2);  end

%% ========================================================================
%  simulate infection transmission on simplicial complex of dim = 3
%  generate NR relaizations
%  ========================================================================
%  pairwise interaction contributes to infection if susceptible node i is 
%  exposed, i.e. if any of its edge-neighbours is infected.
%  A triangle is 'fully active' for node i if other triangle members are 
%  simultaneously infected
%  A tetrahedron is 'fully active' for node i if the three other members
%  are simultaneously infected 

NR = 25; 

I_av = zeros(N,NR);
R_av = zeros(N,NR);
S_av = zeros(N,NR);

tic
parfor run = 1:NR

    I = zeros(n,N);                         % infection state matrix
    R = zeros(n,N);                         % recovery state matrix
    
    I(seed_nodes,1) = 1;                    % set all seed node to 1 initially (infected)

    for t = 1:N-1

        I_new = I(:,t);                     % initiallize next state using current state
        R_new = R(:,t);
        
        xi_t  = randn;                                  
        beta1 = bet_1 + sgma1*xi_t;         % perterbued infection rates
        beta2 = bet_2 + sgma2*xi_t;
        beta3 = bet_3 + sgma3*xi_t;
        beta4 = bet_4 + sgma4*xi_t;
        beta5 = bet_5 + sgma5*xi_t;
        beta6 = bet_6 + sgma6*xi_t;


        for i = 1:n
            
            if I(i,t) == 0                  % susceptible node
        
                % STEP 1: Edge-based (pairwise) contagion transmission
                %   For each edge containing node i, check if the neighbour
                %   is infected. Each infected edge-neighbour contributes β₁
                %   to the infection rate λᵢ.
    
                lambda_i = 0;

                for ix = edge_find{i}'
                    pair = edges(ix,:);
                    nhbr = pair(pair ~= i);
                    if I(nhbr,t) == 1 
                        lambda_i = lambda_i + beta1;
                     end
                end

                % STEP 2: Triangle-based infection
                %   For triangular face containing node i, check if the two
                %   neighbours are infected, the triangle is fully active 
                %   for node i and contribute at rate beta2.

                nFace = 0;

                for ix = tri_find{i}'
                    tri    = triangles(ix,:); 
                    others = tri(tri ~= i);
                    if sum(I(others,t)) == 2 
                        nFace = nFace + 1;
                    end
                end
                lambda_i = lambda_i + beta2*nFace;       % update infection probability
  
                % STEP 3: Tetrahedron-based infection
                %   The tetrahedron is fully active for node i if three
                %   neighbours are infected. 

                nFace = 0;

                for ix = tetra_find{i}'
                    tet    = tetrahedron(ix,:);
                    others = tet(tet ~= i);
                    if sum(I(others, t)) == 3 
                        nFace = nFace + 1;
                    end
                end                                             
                lambda_i = lambda_i + beta3*nFace;     % update infection probability

                % STEP 4: 5-face-based infection
                %   The 5-face is fully active for node i if three
                %   neighbours are infected. 

                nFace = 0;

                for ix = penta_find{i}'
                    pen    = five_face(ix,:);
                    others = pen(pen ~= i);
                    if sum(I(others, t)) == 4 
                        nFace = nFace + 1;
                    end
                end                                             
                lambda_i = lambda_i + beta4*nFace;     % update infection probability

                % STEP 5: 6-face-based infection
                %   The 6-face is fully active for node i if four
                %   neighbours are infected. 

                nFace = 0;

                for ix = hexa_find{i}'
                    hex    = six_face(ix,:);
                    others = hex(hex ~= i);
                    if sum(I(others, t)) == 5 
                        nFace = nFace + 1;
                    end
                end                                             
                lambda_i = lambda_i + beta5*nFace;     % update infection probability

                % STEP 6: 7-face-based infection
                %   The 7-face is fully active for node i if six
                %   neighbours are infected. 

                nFace = 0;

                for ix = hepta_find{i}'
                    hep    = seven_face(ix,:);
                    others = hep(hep ~= i);
                    if sum(I(others, t)) == 6 
                        nFace = nFace + 1;
                    end
                end                                             
                lambda_i = lambda_i + beta6*nFace;     % update infection probability
                
                p_infect = 1 - exp(-lambda_i*dt);
                if rand < p_infect
                    I_new(i) = 1;
                    R_new(i) = 0;
                end

            elseif I(i,t) == 1                  % If the node is infected
                
                % STEP 7: Recovery
                if rand < gamma*dt
                    I_new(i) = 0;
                    R_new(i) = 1;
                end
            end
        end

        % STEP 8: Update states
        I(:,t+1) = I_new;
        R(:,t+1) = R_new;
    end

    I_av(:,run) = mean(I,1);                    % find ensemble infection
    R_av(:,run) = mean(R,1); 
    S_av(:,run) = 1 - I_av(:,run) - R_av(:,run);
end
toc

%% =======================================================================
% Plot sample paths
%=========================================================================
% susceptible 
f = figure;%('visible', 'off');
plot((1:N)*dt,S_av); axis([0 T 0 1]); 
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('$t$','Interpreter','latex', 'FontSize',15);
title('$S(t)$', 'Interpreter','latex', 'FontSize',15);
grid on;
% exportgraphics(f, strcat('S_',dname,'_',fname,'.png'));


% infected 
f = figure;%('visible', 'off');
plot((1:N)*dt,I_av); axis([0 T 0 1]);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('$t$','Interpreter','latex', 'FontSize',15);
title('$I(t)$', 'Interpreter','latex', 'FontSize',15);
grid on;
% exportgraphics(f, strcat('I_',dname,'_',fname,'.png'));


% removed 
f = figure;%('visible', 'off');
plot((1:N)*dt,R_av); axis([0 T 0 1]);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('$t$','Interpreter','latex', 'FontSize',15);
title('$R(t)$', 'Interpreter','latex', 'FontSize',15);
grid on;
% exportgraphics(f, strcat('R_',dname,'_',fname,'.png'));
