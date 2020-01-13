function plot_reference_and_grad_dot_images(im1, im2, figure_handle)
    if nargin < 3
        figure    
    else
        figure(figure_handle)
    end
    
    rmax = max(size(im1, 1), size(im2, 1));
    cmax = max(size(im1, 2), size(im2, 2));
    
    
    subplot(1, 2, 1)
    imagesc(im1)
    colormap(flipud(gray))
    annotate_image(gcf, gca)
    axis([1 cmax 1 rmax])
    set(gca, 'ydir', 'normal')
    title('Reference')
    
    subplot(1, 2, 2)
    imagesc(im2)
    colormap(flipud(gray))
    annotate_image(gcf, gca)
    axis([1 cmax 1 rmax])
    set(gca, 'ydir', 'normal')
    title('Gradient')
    
end