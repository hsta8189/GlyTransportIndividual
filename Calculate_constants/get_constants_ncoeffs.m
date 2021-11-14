function nodes = get_constants_ncoeffs(nodes, nodes_to_keep)
% Reduces the matrix of rate constants nodes to just the nodes in the array
% nodes_to_keep
num_nodes = length(nodes);
if nargin < 2
    nodes_to_keep = [1, num_nodes];
end
mask = ones(num_nodes) - eye(num_nodes); % will use this to get rid of diagonal matrix entries later

n = 1;
num_eliminated = 0;
goal = num_nodes - length(nodes_to_keep);

while num_eliminated < goal
    nodes_old = zeros(size(nodes));
    n = n + 1;
    if n > num_nodes
        n = 2;
    end
    while any(nodes_old ~= nodes, 'all') % sanity check so it wont loop forever
        nodes_old = nodes; 
        for i = 1:num_nodes
            node = nodes(i,:);
            if nnz(node) == n  && ~ismember(i,nodes_to_keep) % we are of the form --node j -- node i -- node k
                syms x [1, num_nodes]
                
                xians =  solve( dot(nodes(:,i), x) - sum(node * x(i))  == 0, x(i)); % concentration of x(i) assuming it is in equilibrium (dxi/dt = 0)
                % get the coeffs of each variable
                syms temp; % dummy variable so that we keep the structure of the terms that dont show up on xians
                xsolns = coeffs(xians + temp * sum(x(:)), fliplr(x));

                xi_content = double(subs(xsolns,temp,0)); % get rid of temporary variable, we now have the coeffs of the states that make up xi

                %neighbouring_nodes = find(xi_content);
                
                % add the number/probability of transporters that would be in xi back into
                % the system while obeying detailed balance
                nodes = nodes + (node .* xi_content') .* mask;
                 nodes(i,:) = 0; 
                 nodes(:,i) = 0;
                 num_eliminated = num_eliminated + 1;
                 
                 fprintf("eliminated %d connections from node %d/%d\n", n,num_eliminated, goal);
                 clear x*;
            end
        end

    end
end


end

