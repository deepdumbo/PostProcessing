function [ output_args ] = saveFigure( handle_keep, fileName )
% saveFigure 
%   Saves figure specified by `handle` as `fileName` in fullscreen
%   as to get around the stupid behavior.
    
    handle = handle_keep;
    set(handle, 'Visible','off');
    screen_size = get(0, 'ScreenSize');
   % origSize = get(handle, 'Position'); % grab original on screen size
    set(handle, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to scren size
    set(handle,'PaperPositionMode','auto') %set paper pos for printing
    saveas(handle, fileName, 'jpg') % save figure
   % set(handle,'Position', origSize) %set back to original dimensions
end