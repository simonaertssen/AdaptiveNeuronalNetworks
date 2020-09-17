function h = phasespaceplot()
    th = 0:pi/50:2*pi;
    unitcircle = [cos(th); sin(th)];
    drawcircle = round(unitcircle);
    
    % Fill the outside
    fill([drawcircle(1,:) flip(unitcircle(1,:))], [drawcircle(2,:) flip(unitcircle(2,:))], 'w')

    % Imaginary unit circle:
    th = 0:pi/50:2*pi;
    h = plot(cos(th), sin(th), '-k', 'LineWidth', 2);
    
    % Plot options:
    grid on; box on; axis square;
end

