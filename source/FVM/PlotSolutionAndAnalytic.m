function PlotSolutionAndAnalytic(nodesX, tout, yout, analyticSol, rows, columns, ...
    extraTitle)

figure;

j = 1;
for i = 3:length(tout)
    
    subplot(2, 2, j);
    j = j + 1;
    
    plot(nodesX, diag(flipud(reshape(yout(:, i), rows, columns) )), 'LineWidth', 2);
    hold all;
    plot(nodesX, diag(flipud(reshape(analyticSol(:, i), rows, columns))), '--r', 'LineWidth', 2);
    
    plotTitle = ['t = ' num2str(tout(i)) ' (' extraTitle ')'];
    title(plotTitle);
    xlabel('x, y');
    ylabel('Solution');

    ylim([0 0.4]);
    
    legend('Numeric', 'Analytic');

    set(findall(gcf,'type','text'), 'fontSize', 12);
        set(gca, 'fontSize', 11);
end
end