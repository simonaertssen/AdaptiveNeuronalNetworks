function fighandle = rasterplot(t, drawthetas, labelfont, ptsze)
    [N, tpts] = size(drawthetas);
    
    [spikesr, spikesc] = ind2sub([N, tpts], find(isnan(drawthetas)));
    fighandle = scatter(t(spikesc), spikesr, ptsze, '.k');
    
    ylim([0, N]); xlim([0, t(end)]);
    ylabel('Neuron index', 'FontSize', labelfont)
end

