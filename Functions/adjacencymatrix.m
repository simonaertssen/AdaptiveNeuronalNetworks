function A = adjacencymatrix(degrees_in, degrees_out)
    N = numel(degrees_in);
    if max(degrees_in) >= N
        warning('Degree too large');
        return
    end
    if nargin == 1
        degrees_out = degrees_in(randperm(N));
    end
    % Test whether the number of elements is the same.
    assert(sum(degrees_in) == sum(degrees_out))
    
    nonzeros = sum(degrees_in);
    xidx = uint16(zeros(nonzeros, 1));
    yidx = uint16(zeros(nonzeros, 1));
    
    choosefrom = uint16(2:N);
    prob_leftout = degrees_out(1);
    probs = degrees_out(choosefrom);

    start = 1;
    for i = 1:N
        num = degrees_in(i);

        % Sample with replacement when enough probabilities are nonzero.
        try
            % sampling without replacement
            chosen = datasample(choosefrom, num, 'Replace', false, 'Weights', probs.^2);
        catch
            % ok, there's not enough probabilities available
            % chosen = maxk(choosefrom, num);
            probs = probs + 1;
            chosen = datasample(choosefrom, num, 'Replace', false, 'Weights', probs+1);
        end

        xidx(start:(start+(num-1))) = i;
        yidx(start:(start+(num-1))) = chosen;

        [~, chosenidx] = ismember(chosen, choosefrom);
        uniqueidx = unique(chosenidx);
        probs(uniqueidx) = probs(uniqueidx) - 1;

        if i ~= N
            choosefrom(i) = i;

            tmp = probs(i);
            probs(i) = prob_leftout;
            prob_leftout = tmp;
        end
        start = start + num;
    end
    
    % Create the matrix as sparse
    A = sparse(xidx, yidx, ones(nonzeros, 1, 'logical'));
    assert(sum(diag(A)) == 0)

    diffrows = degrees_in' - full(sum(A,2))';
    diffcols = degrees_out' - full(sum(A,1));

    N2 = N^2; thresh = 1.0e-6;
    if sum(diffrows)/N2 > thresh
        sprintf(['Adjacency matrix might not be accurate: residue ', num2str(sum(diffrows)/N2)])
    end
    if sum(diffcols)/N2 > thresh
        sprintf(['Adjacency matrix might not be accurate: residue ', num2str(sum(diffcols)/N2)])
    end
end

