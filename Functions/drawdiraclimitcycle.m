function handle = drawdiraclimitcycle()
%     fileID = fopen('diraclimitcycle.bin');
%     z = fread(fileID)
%     fclose(fileID);
    z = readmatrix('diraclimitcycle.dat')
    handle = plot(z(:,1), z(:,2), ':k', 'LineWidth', 1.5);
end

