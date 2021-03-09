function pan = zoomin(ax,areaToMagnify,panPosition)
% AX is a handle to the axes to magnify
% AREATOMAGNIFY is the area to magnify, given by a 4-element vector that 
%      defines the lower-left and upper-right corners of a rectangle [x1 y1 x2 y2]
% PANPOSTION is the position of the magnifying pan in the figure, defined by
%        the normalized units of the figure [x y w h]
%
% Example:
%       x = -5:0.1:5;
%       subplot(3,3,[4 5 7 8])
%       plot(x,cos(x-2),x,sin(x),x,-x-0.5,x,0.1.*x+0.1)
%       ax = gca;
%       area = [-1 -1 1 1];
%       inlarge = subplot(3,3,3);
%       panpos = inlarge.Position;
%       delete(inlarge);
%       inlarge = zoomin(ax,area,panpos);
%       title(inlarge,'Zoom in')
%
% Copyright Eyal Ben-Hur (eyal.ben-hur@mail.huji.ac.il)
%               The Hebrew University in Jerusalem
fig = ax.Parent;
pan = copyobj(ax,fig);
pan.Position = panPosition;
pan.XLim = areaToMagnify([1 3]);
pan.YLim = areaToMagnify([2 4]);
pan.XTick = [];
pan.YTick = [];
rectangle(ax,'Position',...
    [areaToMagnify(1:2) areaToMagnify(3:4)-areaToMagnify(1:2)])
xy = ax2norm(ax,areaToMagnify([1 4;3 2]));
annotation(fig,'line',[xy(1,1) panPosition(1)],...
    [xy(1,2) panPosition(2)+panPosition(4)],'Color','k')
annotation(fig,'line',[xy(2,1) panPosition(1)+panPosition(3)],...
    [xy(2,2) panPosition(2)],'Color','k')
end