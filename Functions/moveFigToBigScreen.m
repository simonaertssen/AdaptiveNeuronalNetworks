function fighandle = moveFigToBigScreen(fighandle)
if nargin == 0
    fighandle = gcf; 
end
MP=get(0,'MonitorPositions');
if size(MP,1)>1
    pos=get(fighandle,'Position');
    pause(0.01); % this seems sometimes necessary on a Mac
    set(fighandle,'Position',[pos(1,2)+MP(2,1:2) pos(3:4)]);
end
end

