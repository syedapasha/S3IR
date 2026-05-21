clc;
% Load the dataset
data = importdata('SFHH.dat'); 


nodes = unique(data(:,2:3));

% Map node values to indices 
nodeMap = containers.Map(nodes, 1:numel(nodes));

edges = [];
for i = 1:size(data, 1)
    edges = [edges; nodeMap(data(i, 2)), nodeMap(data(i, 3))];
end

sortedEdges = sort(edges, 2);  % Sort the nodes within each edge

uniqueEdges = unique(sortedEdges, 'rows');  % Remove duplicate edges

G = graph(uniqueEdges(:,1), uniqueEdges(:,2));
%disp('Filtered unique edges:');
disp(uniqueEdges);
writematrix(uniqueEdges, 'filtered_sfhh.csv');