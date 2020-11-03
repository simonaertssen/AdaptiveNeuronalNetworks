function handle = drawfixeddegreelimitcycle()
    z = readmatrix('diraclimitcycle.dat');
    handle = plot(z(:,1), z(:,2), ':k', 'LineWidth', 1.5);
end

