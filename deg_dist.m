%%========================================================================
%     DEGREE DISTRIBUTION & SIMPLICIAL COMPLEX 
%  =======================================================================
%  Description : Builds a graph (1-simplex) from an edge list, extracts
%                adjacency structure, and enumerates all K-simplices
%                where k=1,2,...,6 present in the dataset. Computes degrees
%                and average node to facet degree.
%                Constructs degree distributions for simplicial complexes
%                (triangles to heptahedra) extracted from a real-world
%                dataset. 
%
%  Input       : edges  — Ex2 matrix of node index pairs [nodeA, nodeB]
%                id     — Nx1 vector of original node IDs
%  Output      : Average Degree
%                Average node-t=facet-degree
%                CSV files representing k-simplices for k=1,2,...,6
%                Histogram describing degree distribution

clc; close all; clear all;

%% ===========================================================
% SECTION 1: Import data
%=============================================================
% The CSV data files describeing two datasets'tij\_InVS' and 'SFHH',
% containg edges sorted without duplicates respectively:
%    'filtered_tijInvs.csv'      each row is a pair   [i, j]
%    'filtered_sfhh.csv'         each row is a pair   [i, j]
%=============================================================

dataset = 'invs';                                   % load datasets
if dataset=='invs'
    data = importdata('filtered_tijInvs.csv');      % sorted, without duplicates 
else dataset=='sfhh'
    data = importdata('filtered_sfhh.csv');         % sorted, without duplicates 
end

id = unique(data);                                  % node IDs

edges = data;


%% =======================================================================
%  SECTION 2 — GRAPH CONSTRUCTION (1-Simplex)
%=========================================================================
%  Build an undirected graph from the edge list and derive the
%  adjacency matrix. The graph represents the 1-skeleton (edges only)
%  on which higher-order simplices will be detected.
% ------------------------------------------------------------------------

G = graph(edges(:,1), edges(:,2));                  % undirected graph
N = numnodes(G);                                    % Total number of nodes
AM = adjacency(G);                                  % adjacency matrix
%  degree(i) = number of edges incident to node i (graph-theoretic degree)
degrees = sum(AM,2);            
avg_degree = mean(degrees);          
disp(['Average node degree: ', num2str(avg_degree)]);

%% =======================================================================
%  SECTION 3 — NEIGHBOUR LIST CONSTRUCTION
%=========================================================================
%  Convert the adjacency matrix into a cell array of neighbour index lists.
%  links{i} contains the column indices of all non-zero entries in row i.
% ------------------------------------------------------------------------

links = cell(N,1);

for i = 1:N
    links{i} = find(AM(i,:));
end

%% =====================================================================
%  SECTION 4 — Construct 2-SIMPLEX  (Triangles)
%=======================================================================
%  Enumerate all closed triangles (3-cliques) in the graph by checking,
%  for every pair of neighbours (j,k) of node i, whether an edge (j,k)
%  also exists. 
%    triangles        — raw order  [i, ni(j), ni(k)]  (insertion order)
%    triangles_sorted — sorted     [0 0 0]     (duplicate guard)
% Compute the average node-to-facet-degree
% ------------------------------------------------------------------------

triangles = [];
triangles_sorted = [0 0 0];


for i = 1:N
    ni = links{i};                        % neighbours of node i
    for j = 1:length(ni)
        for k = j+1:length(ni)
                                          %  Triangle (i, ni(j), ni(k)) is valid iff ni(j)—ni(k) exists
            if AM(ni(j),ni(k))==1 

                candidate = sort([id(i) ni(j) ni(k)]);
                idx = find(triangles_sorted==candidate);
                
                if isempty(idx)           % only add if any permutation is not in the list  
                    triangles = [ triangles; id(i) ni(j) ni(k) ];
                    triangles_sorted = [ triangles_sorted; candidate ];
                end

            end 

        end
    end

end

num_triangles = size(triangles, 1);

% Find nodes forming triangles
idt = unique(triangles);
M = numel(idt);

% Initialize node-to-facet degree for triangles
node_to_facet_degree_tri = zeros(M,1);

% Count occurrences of each node in a triangle
for i = 1:M
    node_to_facet_degree_tri(i) = sum(sum(triangles == idt(i)));
end

% Compute the average node-to-facet degree
avg_facet_degree = sum(node_to_facet_degree_tri)/N;     % divide by M or N?

disp(['Average node-to-facet degree (triangles): ', num2str(avg_facet_degree)]);

disp(triangles);
writematrix(triangles,strcat('triangles_',dataset,'.csv'));
%% ========================================================================
%  SECTION 5 —  Construct 3-SIMPLEX (Tetrahedra)
%==========================================================================
%  Enumerate all closed tetrahedra (4-cliques) in the graph by checking,
%  for every neighbours (j,k,el) of node i, whether edges (j,k) (j,el) (k,el)
%  also exists. 
% Compute the average node-to-facet-degree(tetrahedra)
% ------------------------------------------------------------------------

tetrahedra = [];
tetrahedra_sorted = [0 0 0 0];


for i = 1:N

    ni = links{i};

    for j = 1:length(ni)
        for k = j+1:length(ni)
            for el = k+1:length(ni)

                if AM(ni(j),ni(k))==1 && AM(ni(k),ni(el))==1 && AM(ni(j),ni(el))==1
                                        
                    candidate = sort([id(i) ni(j) ni(k) ni(el)]);
                    idx = find(tetrahedra_sorted==candidate);
                    
                    if isempty(idx)     % only add if any permutation is not in the list  
                        tetrahedra = [ tetrahedra; id(i) ni(j) ni(k) ni(el) ];
                        tetrahedra_sorted = [ tetrahedra_sorted; candidate ];
                    end
                    
                end

            end
        end
    end

end

num_tetrahedrons = size(tetrahedra,1);

% Find nodes forming tetrahedra
idt = unique(tetrahedra);
L = numel(idt);

% Initialize node-to-facet degree (size based on the maximum node index)
node_to_facet_degree = zeros(L,1);

% Count the occurrences of each node in a tetrahedra
for j = 1:L
    node_to_facet_degree(j) = sum(sum(tetrahedra == idt(j)));
end

avg_facet_degree = sum(node_to_facet_degree)/N;     % divide by L or N?

disp(['Average node-to-facet degree (tetrahedra): ', num2str(avg_facet_degree)]);
writematrix(tetrahedra,strcat('tetrahedra_',dataset,'.csv'));
%% ========================================================================
%  SECTION 6 —  Construct 4-SIMPLEX (Pentachora)
%==========================================================================
%  Enumerate all closed Pentachora (5-cliques) in the graph by checking,
%  for every neighbours (j,k,el,m) of node i, whether edges (j,k) (j,el)
%  (k,el)(j,m)(k,m)(el,m) also exists. 
% Compute the average node-to-facet-degree(penta)
% ------------------------------------------------------------------------

pentachora = [];
pentachora_sorted = [0 0 0 0 0];  

for i = 1:N

    ni = links{i};

    for j = 1:length(ni)
        for k = j+1:length(ni)
            for el = k+1:length(ni)
              for m = el+1:length(ni)
                if AM(ni(j),ni(k))==1 && AM(ni(k),ni(el))==1 && AM(ni(j),ni(el))==1 && AM(ni(j),ni(m))==1  && ...
                        AM(ni(m),ni(k))==1  && AM(ni(el),ni(m))==1            
                    candidate = sort([id(i) ni(j) ni(k) ni(el) ni(m)]);
                    idx = find(pentachora_sorted==candidate);
                    
                    if isempty(idx)     % only add if any permutation is not in the list  
                        pentachora = [ pentachora; id(i) ni(j) ni(k) ni(el) ni(m)];
                        pentachora_sorted = [ pentachora_sorted; candidate ];
                    end
                end 
              end

            end
        end
    end

end
num_pentachora = size(pentachora,1);

% Find nodes forming pentachora
idt = unique(pentachora);
R = numel(idt);

% Initialize node-to-facet degree (size based on the maximum node index)
node_to_facet_degree_penta = zeros(R,1);

% Count the occurrences of each node in a pentachora
for r = 1:R
    node_to_facet_degree_penta(r) = sum(sum(pentachora == idt(r)));
end

avg_facet_degree = sum(node_to_facet_degree_penta)/N;     % divide by R or N?

disp(['Average node-to-facet degree (pentachora): ', num2str(avg_facet_degree)]);
writematrix(pentachora,strcat('penta_',dataset,'.csv'));
%% ========================================================================
%  SECTION 7 —  Construct 5-SIMPLEX (hexahedra)
%==========================================================================
%  Enumerate all closed hexahedra(6-cliques) in the graph by checking,
%  for every neighbours (j,k,el,m,p) of node i, whether edges (j,k) (j,el)
%  (k,el)(j,m)(k,m)(el,m) (j,p) (k,p) (el,p) (m,p) also exists. 
%  and compute the average node-to-facet-degree (hexahedra)
% ------------------------------------------------------------------------
hexahedra = [];
hexahedra_sorted = [0 0 0 0 0 0];  

for i = 1:N

    ni = links{i};

    for j = 1:length(ni)
        for k = j+1:length(ni)
            for el = k+1:length(ni)
              for m = el+1:length(ni)
                  for p=m+1:length(ni)

               if AM(ni(j),ni(k))==1 && AM(ni(k),ni(el))==1 && AM(ni(j),ni(el))==1 && ...
                       AM(ni(j),ni(m))==1  && AM(ni(m),ni(k))==1  && AM(ni(el),ni(m))==1 && ...
                       AM(ni(j), ni(p))==1 && AM(ni(k), ni(p))==1 && AM(ni(el),ni(p))==1 && AM(ni(m),ni(p))==1
                    candidate = sort([id(i) ni(j) ni(k) ni(el) ni(m) ni(p)]);
                    idx = find(hexahedra_sorted==candidate);
                    
                    if isempty(idx)     % only add if any permutation is not in the list  
                        hexahedra = [hexahedra; id(i) ni(j) ni(k) ni(el) ni(m) ni(p)];
                        hexahedra_sorted = [ hexahedra_sorted; candidate ];
                    end
               end
                 end 
              end

            end
        end
    end

end
num_hexahedra = size(hexahedra,1);

% Find nodes forming hexahedra
idt = unique(hexahedra);
S = numel(idt);

% Initialize node-to-facet degree (size based on the maximum node index)
node_to_facet_degree_hexa = zeros(S,1);

% Count the occurrences of each node in a hexahedra
for s = 1:S
    node_to_facet_degree_hexa(s) = sum(sum(hexahedra == idt(s)));
end

avg_facet_degree = sum(node_to_facet_degree_hexa)/N;     % divide by S or N?

disp(['Average node-to-facet degree (hexahedra): ', num2str(avg_facet_degree)]);
%disp(hexahedra);
writematrix(hexahedra,strcat('hexa_',dataset,'.csv'));


%% ========================================================================
%  SECTION 8 —  Construct 6-SIMPLEX (heptahedra)
%==========================================================================
%  Enumerate all closed heptahedra(7-cliques) in the graph by checking,
%  for every neighbours (j,k,el,m,p,q) of node i, whether edges (j,k) (j,el)
%  (k,el)(j,m)(k,m)(el,m) (j,p) (k,p) (el,p) (m,p) (j,q) (k,q) (el,q) (m,q)
% (p,q) also exists. 
%  and compute the average node-to-facet-degree(heptahedra)
% ------------------------------------------------------------------------
heptahedra = [];
heptahedra_sorted = [0 0 0 0 0 0 0];  

for i = 1:N

    ni = links{i};

    for j = 1:length(ni)
        for k = j+1:length(ni)
            for el = k+1:length(ni)
              for m = el+1:length(ni)
                  for p=m+1:length(ni)
                      for q=p+1:length(ni)

               if AM(ni(j),ni(k))==1 && AM(ni(k),ni(el))==1 && AM(ni(j),ni(el))==1 && ...
                       AM(ni(j),ni(m))==1  && AM(ni(m),ni(k))==1  && AM(ni(el),ni(m))==1 && ...
                       AM(ni(j), ni(p))==1 && AM(ni(k), ni(p))==1 && AM(ni(el),ni(p))==1 && AM(ni(m),ni(p))==1 &&...
                       AM(ni(j), ni(q))==1 && AM(ni(k), ni(q))==1 && AM(ni(el), ni(q))==1 && ...
                       AM(ni(m), ni(q))==1 && AM(ni(p), ni(q))==1 
                    candidate = sort([id(i) ni(j) ni(k) ni(el) ni(m) ni(p) ni(q)]);
                    idx = find(heptahedra_sorted==candidate);
                    
                    if isempty(idx)     % only add if any permutation is not in the list  
                        heptahedra = [heptahedra; id(i) ni(j) ni(k) ni(el) ni(m) ni(p) ni(q)];
                        heptahedra_sorted = [ heptahedra_sorted; candidate ];
                    end
               end
                      end
                 end 
              end

            end
        end
    end

end
num_heptahedra = size(heptahedra,1);

% Find nodes forming heptahedra
idt = unique(heptahedra);
S = numel(idt);

% Initialize node-to-facet degree (size based on the maximum node index)
node_to_facet_degree_hepta = zeros(S,1);

% Count the occurrences of each node in a heptahedra
for s = 1:S
    node_to_facet_degree_hepta(s) = sum(sum(heptahedra == idt(s)));
end

avg_facet_degree = sum(node_to_facet_degree_hepta)/N;     % divide by S or N?

disp(['Average node-to-facet degree (heptahedra): ', num2str(avg_facet_degree)]);
writematrix(heptahedra,strcat('hepta_',dataset,'.csv'));

%% ========================================================================
%  SECTION 9 —  Construct Histogram
%==========================================================================
%  Plot the degree distribution histogram for each simplex
% ------------------------------------------------------------------------
nbins = 50;

% plot degree distribution 
f = figure('visible', 'off');
histogram(degrees, nbins, 'Normalization','probability'); 
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('degree','Interpreter','latex', 'FontSize',15);
title('Histogram', 'Interpreter','latex', 'FontSize',15);
grid on;
exportgraphics(f, strcat('deg_',dataset,'.png'));


nbins = 5;
% plot node-to-face (triangle) degree distribution 
f = figure('visible', 'off');
histogram(node_to_facet_degree_tri, nbins, 'Normalization','probability'); 
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('Node to 3-face degree','Interpreter','latex', 'FontSize',15);
title('Histogram', 'Interpreter','latex', 'FontSize',15);
grid on;
exportgraphics(f, strcat('TRdeg_',dataset,'.png'));


nbins = 4;
% plot node-to-face (tetrhedra) degree distribution 
f = figure('visible', 'off');
histogram(node_to_facet_degree, nbins, 'Normalization','probability'); 
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('Node to 4-face degree','Interpreter','latex', 'FontSize',15);
title('Histogram', 'Interpreter','latex', 'FontSize',15);
grid on;
exportgraphics(f, strcat('TTdeg_',dataset,'.png'));

nbins = 5;
% plot node-to-face (pentachora) degree distribution 
f = figure('visible', 'off');
histogram(node_to_facet_degree_penta, nbins, 'Normalization','probability'); 
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('Node to 5-face degree','Interpreter','latex', 'FontSize',15);
title('Histogram', 'Interpreter','latex', 'FontSize',15);
grid on;
exportgraphics(f, strcat('PENTdeg_',dataset,'.png'));

nbins = 5;
% plot node-to-face (hexahedra) degree distribution 
f = figure('visible', 'off');
histogram(node_to_facet_degree_hexa, nbins, 'Normalization','probability'); 
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('Node to 6-face degree','Interpreter','latex', 'FontSize',15);
title('Histogram', 'Interpreter','latex', 'FontSize',15);
grid on;
exportgraphics(f, strcat('HEXdeg_',dataset,'.png'));


nbins = 5;
%plot node-to-face (Heptahedra) degree distribution 
f = figure('visible', 'off');
histogram(node_to_facet_degree_hepta, nbins, 'Normalization','probability'); 
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('Node to 7-face degree','Interpreter','latex', 'FontSize',15);
title('Histogram', 'Interpreter','latex', 'FontSize',15);
grid on;
exportgraphics(f, strcat('HEPdeg_',dataset,'.png'));