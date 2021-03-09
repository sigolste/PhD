function normxy = ax2norm(ax,xy)
% This function converts the axis unites to the figure normalized unites
%   AX is a handle to the figure
%   XY is a n-by-2 matrix, where the first column is the x values and the
%       second is the y values
%   ANXY is a matrix in the same size of XY, but with all the values
%       converted to normalized units
%
% Copyright Eyal Ben-Hur (eyal.ben-hur@mail.huji.ac.il)
%               The Hebrew University in Jerusalem
pos = ax.Position;
%  white area * ((value - axis min)  / axis length)  + gray area
normx = pos(3)*((xy(:,1)-ax.XLim(1))./range(ax.XLim))+ pos(1);
normy = pos(4)*((xy(:,2)-ax.YLim(1))./range(ax.YLim))+ pos(2);
normxy = [normx normy];
end