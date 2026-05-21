clc
rng(100,"twister")

%% Import data
edges       = readmatrix('edges_sfhh.csv');
triangles   = readmatrix('triangles_sfhh.csv');
tetrahedron = readmatrix('tetrahedra_sfhh.csv');
n = max([edges(:); triangles(:); tetrahedron(:)]);

T  = 10;
dt = 0.001;
N  = T/dt;


%% Pre-compute neighbour array 
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

%% initial infected seed

seed_nodes = randperm(n, 2);             % draw seed nodes


NR=1; 

f = figure('visible', 'on');
 hold on;
for run=1:NR 
%%  Edge-based infection (D=1)
seed_nodes = randperm(n, 2);
gamma  = 0.5;
    beta_1 = 0.1;   % outbreak will die out
    sigma1 = 0.01;  

I1 = zeros(n,N);
R1 = zeros(n,N);
I1(seed_nodes, 1) = 1;
R1(:,1) = 0;

for t = 1:N-1
    I1_new = I1(:, t);
    R1_new = R1(:, t);
    xi_t   = randn;
    beta1  = beta_1 + sigma1*xi_t;         % perturbed, 
%     beta2  = beta_2 + sigma2*xi_t;
%     beta3  = beta_3 + sigma3*xi_t;

    for i = 1:n
        if I1(i,t) == 0                             % susceptible node
            lambda_i = 0;
            % STEP1: Edge-based infection
            for e = edge_find{i}'
                pair = edges(e, :);
                nhbr = pair(pair ~= i);
                if I1(nhbr, t) == 1
                    lambda_i = lambda_i + beta1;
                end
            end
            p_infect = 1 - exp(-lambda_i*dt);
            if rand < p_infect
                I1_new(i) = 1;
                R1_new(i) = 0;
            end

        elseif I1(i,t) == 1
            % STEP4: Recovery
            if rand < gamma*dt
                I1_new(i) = 0;
                R1_new(i) = 1;
            end
        end
    end
    % STEP5: Update states
    I1(:, t+1) = I1_new;
    R1(:, t+1) = R1_new;
end

I1_t = mean(I1, 1);                               % Find average infection
R1_t = mean(R1, 1);
S1_t = 1 - I1_t - R1_t;

%% infection in 2-simplex (D=2)
seed_nodes = randperm(n, 2);
gamma  = 0.4;
    beta_1 = 0.1; beta_2 = 0.2;   % outbreak will die out
    sigma1 = 0.01;  sigma2 = 0.02;  

I2 = zeros(n,N);
R2 = zeros(n,N);
I2(seed_nodes, 1) = 1;
R2(:,1) = 0;

for t = 1:N-1
    I2_new = I2(:, t);
    R2_new = R2(:, t);
    xi_t   = randn;
    beta1  =  beta_1 + sigma1*xi_t;
    beta2  =  beta_2 + sigma2*xi_t;
    

    for i = 1:n
        if I2(i,t) == 0                             % susceptible node
            lambda_i = 0;
            % STEP1: Edge-based infection
            for e = edge_find{i}'
                pair = edges(e, :);
                nhbr = pair(pair ~= i);
                if I2(nhbr, t) == 1 
                    lambda_i = lambda_i + beta1;
                end
            end
            % STEP2: Triangle-based infection
            lambda_tri = 0;
            for k = tri_find{i}'
                tri    = triangles(k, :);
                others = tri(tri ~= i);
                if sum(I2(others, t)) == 2 
                    lambda_tri = lambda_tri + 1;
                end
            end
            lambda_i = lambda_i + beta2 * lambda_tri;      
            p_infect = 1 - exp(-lambda_i*dt);
            if rand < p_infect
                I2_new(i) = 1;
                R2_new(i) = 0;
            end

        elseif I2(i,t) == 1
            % STEP4: Recovery
            if rand < gamma*dt
                I2_new(i) = 0;
                R2_new(i) = 1;
            end
        end
    end
    % STEP5: Update states
    I2(:, t+1) = I2_new;
    R2(:, t+1) = R2_new;
end

I2_t = mean(I2, 1);                                  % find average infection
R2_t = mean(R2, 1);
S2_t = 1 - I2_t - R2_t;

%% infection 3-simplex (D=3)
seed_nodes = randperm(n, 3);
    gamma  = 0.3;
    beta_1 = 0.1; beta_2 = 0.2; beta_3 = 0.3;  % outbreak will die out
    sigma1 = 0.01;  sigma2 = 0.02;  sigma3 = 0.03;
    
I3 = zeros(n,N);
R3 = zeros(n,N);
I3(seed_nodes, 1) = 1;
R3(:,1) = 0;

for t = 1:N-1
    I3_new = I3(:, t);
    R3_new = R3(:, t);
    xi_t  = randn;
    beta1 =  beta_1 + sigma1*xi_t;
    beta2 =  beta_2 + sigma2*xi_t;
    beta3 =  beta_3 + sigma3*xi_t;

    for i = 1:n
        if I3(i,t) == 0                              % susceptible node
            lambda_i = 0;
            % STEP1: Edge-based infection
            for e = edge_find{i}'
                pair = edges(e, :);
                nhbr = pair(pair ~= i);
                if I3(nhbr, t) == 1 
                    lambda_i = lambda_i + beta1;
                end
            end
            % STEP2: Triangle-based infection
            lambda_tri = 0;
            for k = tri_find{i}'
                tri    = triangles(k, :);
                others = tri(tri ~= i);
                if sum(I3(others, t)) == 2 
                    lambda_tri = lambda_tri + 1;
                end
            end
            lambda_i = lambda_i + beta2 * lambda_tri;      

            % STEP3: Tetrahedron-based infection
            lambda_tet = 0;
            for h = tetra_find{i}'
                tet    = tetrahedron(h, :);
                others = tet(tet ~= i);
                if sum(I3(others, t)) >= 2 
                    lambda_tet = lambda_tet + 1;
                end
            end
            lambda_i = lambda_i + beta3 * lambda_tet;
            p_infect = 1 - exp(-lambda_i*dt);
            if rand < p_infect
                I3_new(i) = 1;
                R3_new(i) = 0;
            end

        elseif I3(i,t) == 1
            % STEP4: Recovery
            if rand < gamma*dt
                I3_new(i) = 0;
                R3_new(i) = 1;
            end
        end
    end
    % STEP5: Update states
    I3(:, t+1) = I3_new;
    R3(:, t+1) = R3_new;
end

I3_t = mean(I3, 1);                              % Find average infection
R3_t = mean(R3, 1);
S3_t = 1 - I3_t - R3_t;



%% construct plots

% plot d=1 

plot(dt:dt:T, I1_t);

% plot d=2

plot(dt:dt:T, I2_t);

% plot d=3 
 plot(dt:dt:T, I3_t);

ylim([0 1])
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 15);
title('$I(t)$',  'Interpreter', 'latex', 'FontSize', 15);
grid on;
box on;
%exportgraphics(f, strcat(fname,'_invs.png'));
legend('$D$=1','$D$=2','$D$=3','Interpreter','latex','FontSize',10,'location','southeast');
exportgraphics(f, 'comp_sfhh.png');
end
hold off
