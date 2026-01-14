function cycles_counts = count_cycles_simpleDFS(data,sets, nodes)
    % data: N x 6 cell array with values 'x','y','z','w'

    node_map = containers.Map(nodes,1:4);

    Nparticipants = size(data,1);
    cycles_counts = zeros(Nparticipants,1);

    for p = 1:Nparticipants
        % Build graph edges for participant p
        edges_s = zeros(6,1);
        edges_t = zeros(6,1);
        for c = 1:6
            pair = sets{c};
            choice = data{p,c};
            if ~ismember(choice,pair)
                error('Choice "%s" not in pair %s at participant %d, col %d',choice,strjoin(pair,','),p,c);
            end
            other = pair{~strcmp(pair,choice)};
            edges_s(c) = node_map(choice);
            edges_t(c) = node_map(other);
        end

        % Build adjacency list
        adj = cell(4,1);
        for i=1:4
            adj{i} = [];
        end
        for i=1:6
            adj{edges_s(i)} = [adj{edges_s(i)}, edges_t(i)];
        end

        all_cycles = {};

        for startNode=1:4
            visited = false(4,1);
            path = [];
            cycles_found = dfs_cycles(startNode, startNode, adj, visited, path);
            all_cycles = [all_cycles; cycles_found];
        end

        % Remove duplicate cycles (considering rotation and reversal)
        unique_cycles_str = unique(cellfun(@canonical_cycle_str, all_cycles, 'UniformOutput', false));
        cycles_counts(p) = numel(unique_cycles_str);
    end
end

function cycles_found = dfs_cycles(currentNode, startNode, adj, visited, path)
    path = [path currentNode];
    visited(currentNode) = true;
    cycles_found = {};

    for nxt = adj{currentNode}
        if ~visited(nxt)
            cycles_found = [cycles_found; dfs_cycles(nxt, startNode, adj, visited, path)];
        elseif nxt == startNode && numel(path) > 1
            % Found a cycle!
            cycles_found = [cycles_found; {path}];
        end
    end
end

function cstr = canonical_cycle_str(c)
    % Normalize cycle c to remove duplicates caused by rotations or reversal

    n = numel(c);
    % Remove duplicated end nodes if present
    if c(1) == c(end)
        c(end) = [];
        n = n - 1;
    end

    % All rotations:
    c_double = [c c];
    candidates = cell(1,n);
    for i=1:n
        candidates{i} = c_double(i:i+n-1);
    end

    % Minimal rotation:
    min_seq = candidates{1};
    for i=2:n
        if lexless(candidates{i}, min_seq)
            min_seq = candidates{i};
        end
    end
    % Check reversed minimal rotation:
    crev = flip(min_seq);
    if lexless(crev, min_seq)
        min_seq = crev;
    end

    % Convert min_seq to string
    cstr = sprintf('%d-', min_seq);
end

function bool = lexless(a,b)
    % Returns true if vector a lex less than b
    idx = find(a ~= b, 1);
    if isempty(idx)
        bool = false; % equal
    else
        bool = (a(idx) < b(idx));
    end
end