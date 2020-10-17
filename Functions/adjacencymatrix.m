function A = adjacencymatrix(degrees_in, degrees_out)
    N = numel(degrees_in);
    assert(sum(degrees_in) == sum(degrees_out));
    
    if max(degrees_in) >= N
        error('Degree too large');
    end
    if all(degrees_in == N - 1) && all(degrees_in == N - 1)
        A = ones(N) - eye(N);
        disp('Found exact A after 1 try');
        return
    end
    numnonzeros = sum(degrees_in);
    
    % Test for laptop version or other:
    if version('-release') == "2020a"
        numtype = 'uint16';
    else
        numtype = 'double';
    end
    
    tries = 0;
    while tries < 10  
        tries = tries + 1;
        
        xidx = zeros(numnonzeros, 1, numtype);
        yidx = zeros(numnonzeros, 1, numtype);
        idxidx = cumsum([1; degrees_in]); % For indexing the idx and yidx vector

        probs = degrees_out;
        rowpermutation = randperm(N);
        for i = 1:N
            rowindex = rowpermutation(i);
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
            if tries == 1
                disp(['Found exact A after ', num2str(tries), ' try']);
            else
                disp(['Found exact A after ', num2str(tries), ' tries']);
            end
            return
        end
    end
    error(['A might not be accurate: residue of ', num2str(sum(abs(diffrows))/numnonzeros * 100)]);
end

