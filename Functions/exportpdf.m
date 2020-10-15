function fighandle = exportpdf(export, fighandle, figname)
    removewhitspace();
    if export == true
        set(fighandle,'PaperOrientation','landscape');
        set(fighandle,'PaperUnits','normalized');
        set(fighandle,'PaperPosition', [0 0 5 10]);
        print(fighandle, strcat('../Figures/', figname), '-dpdf');
    end
end

