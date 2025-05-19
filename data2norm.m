function norm_coords = data2norm(ax, data_coords)
    % Converts data coordinates to normalized figure coordinates
    ax_units = get(ax, 'Units');
    set(ax, 'Units', 'normalized');
    axpos = get(ax, 'Position');
    set(ax, 'Units', ax_units);

    xlim = get(ax, 'XLim');
    ylim = get(ax, 'YLim');

    norm_x = axpos(1) + (data_coords(1) - xlim(1)) / (xlim(2) - xlim(1)) * axpos(3);
    norm_y = axpos(2) + (data_coords(2) - ylim(1)) / (ylim(2) - ylim(1)) * axpos(4);

    norm_coords = [norm_x, norm_y];
end