function triplot3(n4ed,x,y,z)
% plots a triangular mesh in 2D (red triangles)
%       and its corresponding solution in 3D (blue triangles)
%
% Use: triplot3(n4e(:,[1 2 3 1]),c4n(:,1),c4n(:,2),u)
%  or  triplot3(n4ed,c4n(:,1),c4n(:,2),u)
%  n4e - nodes for elements: each triangle corresponds to one row of three vertices.
%  c4n - coordinates for nodes: each vertex corresponds to one row of coordinates.
%        x=c4n(:,1) consists of the x-coordinates of the vertices
%        y=c4n(:,2) consists of the y-coordinates of the vertices
%  z=u - is the solution, it consists of the discrete solution at each vertex
%
% n4ed - nodes for edges: each edge corresponds to one row of two vertices.


% Copyright 2007 Joscha Gedicke, Hella Rabus
%

hold on; 
for i = 1 : size(n4ed,1)
    curEdge =n4ed(i,:);
    plot3(x(curEdge),y(curEdge),zeros(size(z(curEdge))),'r',"linewidth", 2);
    plot3(x(curEdge),y(curEdge),z(curEdge),'b',"linewidth", 2);
end
hold off

%view(-20,30)
%for d=-20:10:80
%    view ( d, 30 );
%    pause(0.5)
%end
