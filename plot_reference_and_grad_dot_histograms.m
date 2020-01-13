function plot_reference_and_grad_dot_histograms(x1, x2, xlabel_string, figure_handle)
    if nargin < 4
        figure    
    else
        figure(figure_handle)
    end
    
    min_all = min(min(x1), min(x2));
    max_all = max(max(x1), max(x2));
    
    xmin = min_all - (max_all - min_all) * 0.2;
    xmax = max_all + (max_all - min_all) * 0.2;
    
    subplot(1, 2, 1)
    histogram(x1)
    xlim([xmin, xmax])
    xlabel(xlabel_string)
    title('Reference')
    
    subplot(1, 2, 2)
    histogram(x2)
    xlim([xmin, xmax])
    xlabel(xlabel_string)
    title('Gradient')   

    set(gcf, 'Position', [360   632   654   268])
end