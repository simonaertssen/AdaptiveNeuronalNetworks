function A = adjacencymatrix(degrees_in, degrees_out)
    N = numel(degrees_in);
    if max(degrees_in) >= N
        error('Degree too large');
    end
    numnonzeros = sum(degrees_in);
    
    % Test for laptop version or other:
    if version('-release') == "2020a"
        numtype = 'uint16';
    else
        numtype = 'double';
    end
    
    tries = 0;
    while tries < 5  
        tries = tries + 1;
        
        xidx = zeros(numnonzeros, 1, numtype);
        yidx = zeros(numnonzeros, 1, numtype);
        idxidx = cumsum([1; degrees_in]); % For indexing the idx and yidx vector

        probs = degrees_out;

        for i = 1:N
            rowindex = i;
            numelements = degrees_in(rowindex);
            if numelements == 0
                continue
            end
            prob_leftout = probs(rowindex); % Take out diagonal element
            probs(rowindex) = -1;

            % Permutation makes the implementation quite robust:
            % Don't just sample the first maximum elements
            probsperm = randperm(N);
            % [~, probsperminv] = sort(probsperm); % Inverse not necessary?
            [~, chosenidx] = maxk(probs(probsperm), numelements);
            chosenidx = probsperm(chosenidx);

            indices = idxidx(rowindex):idxidx(rowindex+1)-1;
            xidx(indices) = rowindex;
            yidx(indices) = chosenidx;

            probs(chosenidx) = probs(chosenidx) - 1;

            % Reset the probability vector:
            probs(rowindex) = prob_leftout;

        end

        % Create the matrix as sparse
        if gpuDeviceCount > 0
            A = zeros(N, N, 'double');
            A(sub2ind([N, N], xidx, yidx)) = 1;
        else 
            A = sparse(xidx, yidx, ones(numnonzeros, 1, 'logical'));
        end
        A(N,N) = 0; % Make it an N x N matrix
        assert(sum(diag(A)) == 0);

        diffrows = degrees_in' - full(sum(A,2))';
        diffcols = degrees_out' - full(sum(A,1));
        
        if ~any(diffrows) && ~any(diffcols) % Then we have it all right
            disp(['Found exact A after ', num2str(tries), 'tries']);
            return
        end
    end
    
    sprintf('Adjacency matrix might not be accurate: residue of %d', sum(diffrows)/numnonzeros * 100)
end

