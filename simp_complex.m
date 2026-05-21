%==========================================================================
% Construction of Simplicial Complex of dimension D=3 and Compute average 
% node-to-facet-degree.
%==========================================================================
% This script constructs a simplicial complex up to dimension 3 from a real 
% datasets consisting of edge-list and compute avarage node to facet degree
% for each dimension D=1,2,3.
% -------------------------------------------------------------------------
% A simplicial complex is a collection of simplices closed under the
% face relation — every face of a simplex in the complex must also
% belong to the complex.
%--------------------------------------------------------------------------
% Simplex dimensions:
%   1-simplex : edge    (line segment between 2 vertices)
%   2-simplex : triangle (filled triangle among 3 vertices)
%   3-simplex : tetrahedron (filled region among 4 vertices)


clc; close all; clear all;
%% ===========================================================
% SECTION 1: Import data
%=============================================================
% The CSV data files describeing two datasets'tij\_InVS' and 'SFHH',
% containg edges sorted without duplicates respectively:
%    'filtered_tijInvs.csv'      each row is a pair   [i, j]
%    'filtered_sfhh.csv'         each row is a pair   [i, j]
%=============================================================
dataset = 'invs';
%dataset = 'sfhh';
                                                   % import data
if dataset=='invs'
    data = importdata('filtered_tijInvs.csv');      % sorted, without duplicates 
else dataset=='sfhh'
    data = importdata('filtered_sfhh.csv');         % sorted, without duplicates 
end
%% =========================================================================
%SECTION 2: construct (1-simplex) and compute average node degree
%==========================================================================
% A 1-simplex {i, j} is a direct connection between two vertices.
% node IDs representing an undirected edge
% % Average node degree across the network

id = unique(data);                                  % node IDs
edges = data;                                       %  edge list
G = graph(edges(:,1), edges(:,2));                  % undirected graph
N = numnodes(G);                                    % total no of nodes in the graph
AM = adjacency(G);                                  % adjacency matrix
degrees = sum(AM,2);                                % Degree of each node
avg_degree = mean(degrees);                         % Average node degree 

disp(['Average node degree: ', num2str(avg_degree)]);
%% =========================================================================
%SECTION 3:Construct (2-simplex) and compute Average node-to-facet degree (triangles)
%==========================================================================
% A 2-simplex {i, j, k} exists when every node i is connected to two other
% nodes j, k making a triangle with edges i-j, i-k,j-k.
% Equivalently, this finds all triangles in the underlying graph.
links = cell(N,1);

for i = 1:N
    links{i} = find(AM(i,:));        % Indices of nodes connected to node i
end

triangles = [];                     % Each row: [i, j, k] representing a triangle
triangles_sorted = [0 0 0];         % Sort unique triangles

for i = 1:N                         % Loop over all nodes to identify triangles
    ni = links{i};                  % neighbors of node i
    for j = 1:length(ni)
        for k = j+1:length(ni)
                                   % Check if neighbors ni(j) and ni(k) are also connected
                                   % if yes i,j,k forms a triangle
            if AM(ni(j),ni(k))==1 

                candidate = sort([id(i) ni(j) ni(k)]);
                idx = find(triangles_sorted==candidate);
                
                if isempty(idx)     % only add if any permutation is not in the list  
                    triangles = [ triangles; id(i) ni(j) ni(k) ];
                    triangles_sorted = [ triangles_sorted; candidate ];
                end
            end 
        end
    end
end

num_triangles = size(triangles, 1);        % Total number of triangles                               
idt = unique(triangles);                   % Find nodes forming triangles
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

%disp(triangles);

%% =========================================================================
%SECTION 4: construct (3-simplex) and compute Average node-to-facet degree (tetrahedra)
%==========================================================================
% A 3-simplex {i, j, k, l} exists when every node i is connected to three other
% nodes j, k and l making a tetrahedra with edges i-j, i-k, i-l, j-k, j-l, k-l.
% Equivalently, this finds all tetrahedra in the underlying graph.
tetrahedra = [];                    % Each row: [i, j, k, l] representing a tetrahedra       
tetrahedra_sorted = [0 0 0 0];      % Sort unique tetrahedra

for i = 1:N                         % loop over all nodes to identify tetrahedra
    ni = links{i};                  % neighbours of node i
    for j = 1:length(ni)
        for k = j+1:length(ni)
            for el = k+1:length(ni)
                                   % check if all neighbours are connected
                                   % with each other if yes it's a
                                   % tetrahedra
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

num_tetrahedrons = size(tetrahedra,1); % tptal number of tetrahedra

                                       % Find nodes forming tetrahedra
idt = unique(tetrahedra);
L = numel(idt);
                                       % Initialize node-to-facet degree 
                                       % (size based on the maximum node index)
node_to_facet_degree = zeros(L,1);
                                       % Count the occurrences of each node in a tetrahedra
for j = 1:L
    node_to_facet_degree(j) = sum(sum(tetrahedra == idt(j)));
end
                                        % Compute average node to facet degree
avg_facet_degree = sum(node_to_facet_degree)/N; % divide by L or N?

disp(['Average node-to-facet degree (tetrahedra): ', num2str(avg_facet_degree)]);

%% ========================================================================
% SECTION 5: Plot simplices
%==========================================================================
%  Visualise simplicial complex derived from a real
%  datasets across four figures:
%
%   Figure 1 — 1-simplices : the raw graph (vertices + edges only)
%   Figure 2 — 2-simplices : triangular faces shaded yellow
%   Figure 3 — 3-simplices : tetrahedral faces shaded red
%   Figure 4 — Combined    : full simplicial complex, all dimensions
%-----------------------------------------------------------------------
%STEP1: Build graph from edges  
%------------------------------------------------------------------------
G = graph(edges(:,1), edges(:,2));   
                                     % Get layout coordinates (2D)
coords2D = plot(G,'Layout','force'); % force-directed layout
x = coords2D.XData;
y = coords2D.YData;
                                     %  Extend to 3D (add random z)
z = rand(size(x));  
vertices = [x(:), y(:), z(:)];

%-----------------------------------------------------------------------
%STEP2: Show 1-simplex( pairwise connections) 
%------------------------------------------------------------------------

figure;
hold on; 
for e = 1:size(edges,1)
    v = vertices(edges(e,:),:);
    plot3(v(:,1), v(:,2), v(:,3), 'k');
end

scatter3(vertices(:,1), vertices(:,2), vertices(:,3), 40, 'b', 'filled', 'MarkerEdgeColor','k');
axis off;

%-----------------------------------------------------------------------
%STEP3: Show 2-simplex( triangular connections) 
%------------------------------------------------------------------------
figure; hold on; 
for e = 1:size(edges,1)
    v = vertices(edges(e,:),:);
    plot3(v(:,1), v(:,2), v(:,3), 'k');
end

for t = 1:size(triangles,1)
    patch('Faces',1:3,'Vertices',vertices(triangles(t,:),:), ...
          'FaceColor','yellow','FaceAlpha',0.5);
end
scatter3(vertices(:,1), vertices(:,2), vertices(:,3), 40, 'b', 'filled', 'MarkerEdgeColor','k');
axis off;

%-----------------------------------------------------------------------
%STEP4: Show 3-simplex (tetrahedral connections)
%------------------------------------------------------------------------
figure; hold on;
for e = 1:size(edges,1)
    v = vertices(edges(e,:),:);
    plot3(v(:,1), v(:,2), v(:,3), 'k');
end

faces = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
for tet = 1:size(tetrahedra,1)
    v = vertices(tetrahedra(tet,:),:);
    patch('Faces',faces,'Vertices',v, ...
          'FaceColor','red','FaceAlpha',0.1,'EdgeColor','none');
end
scatter3(vertices(:,1), vertices(:,2), vertices(:,3), 40, 'b', 'filled', 'MarkerEdgeColor','k');
axis off;
%-----------------------------------------------------------------------
%STEP5: Plot simplicial complex ( paie+triangles+tetrahedra) 
%------------------------------------------------------------------------
figure; hold on; 

% Edges
for e = 1:size(edges,1)
    v = vertices(edges(e,:),:);
    plot3(v(:,1), v(:,2), v(:,3), 'k');
end

% Triangles
for t = 1:size(triangles,1)
    patch('Faces',1:3,'Vertices',vertices(triangles(t,:),:), ...
          'FaceColor','yellow','FaceAlpha',0.3,'EdgeColor','none');
end

% Tetrahedra
faces = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
for tet = 1:size(tetrahedra,1)
    v = vertices(tetrahedra(tet,:),:);
    patch('Faces',faces,'Vertices',v, ...
          'FaceColor','red','FaceAlpha',0.1,'EdgeColor','none');
end
scatter3(vertices(:,1), vertices(:,2), vertices(:,3), 40, 'b', 'filled', 'MarkerEdgeColor','k');
axis off;
