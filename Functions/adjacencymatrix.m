function A = adjacencymatrix(degrees_in, degrees_out)
    initarray = make_GPUhandle();
    N = numel(degrees_in);
    if max(degrees_in) >= N
        error('Degree too large');
    end
    if nargin == 1
        degrees_out = degrees_in(randperm(N));
    end
    nonzeros = sum(degrees_in);
    % Test for laptop version or other:
    
    if version('-release') == "2020a"
        numtype = 'uint16';
    else
        numtype = 'single';
    end
    xidx = initarray(zeros(nonzeros, 1, numtype));
    yidx = initarray(zeros(nonzeros, 1, numtype));

    choosefrom = cast(2:N,numtype);
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
            chosen = datasample(choosefrom, num, 'Replace', false, 'Weights', probs);
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
    if gpuDeviceCount > 0
        A = initarray(zeros(N, N, 'logical'));
        A(sub2ind([N, N], xidx, yidx)) = 1;
    else 
        A = initarray(sparse(xidx, yidx, ones(nonzeros, 1, 'logical')));
    end
    assert(sum(diag(A)) == 0);

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

