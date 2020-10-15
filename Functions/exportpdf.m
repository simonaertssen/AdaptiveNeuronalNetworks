function fighandle = exportpdf(fighandle, figname, export)
    if export == true
        removewhitspace();
        set(fighandle,'PaperPositionMode','auto', 'PaperOrientation', 'landscape');
        set(gca,'units','centimeters');
        pos = get(gca,'Position');
        ti = get(gca,'TightInset');

        set(fighandle, 'PaperUnits','centimeters');
        set(fighandle, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
        set(fighandle, 'PaperPositionMode', 'manual');
        set(fighandle, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

        print(fighandle, figname, '-dpdf', '-bestfit');
    end
end

