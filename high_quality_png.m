function high_quality_png(filename)
    axis off;
    exportgraphics(gca, filename, 'Resolution', 500);
end