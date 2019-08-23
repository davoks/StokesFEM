function ShowFuncP1 (element, coordinate, u)

% SHOW displays the solution of the finite element computation.
%  Parameters:
%
%    Input, integer ELEMENTS3(N3,3), the nodes that make up each triangle.
%    Input, real COORDINATES(N,1:2), the coordinates of each node.
%    Input, real U(N), the finite element coefficients which represent the solution.
%    There is one coefficient associated with each node.
%
%  Display the information associated with triangular elements.
%
  hold on
%  for j=1:size(element,1)
%    trisurf([1 2 3], coordinate(element(j,:),1), coordinate(element(j,:),2), ...
%	    u(element(j,:)), 'facecolor', 'interp');
%  end
  trisurf(element, coordinate(:,1), coordinate(:,2), u, 'facecolor', 'interp');
  hold off
  view (-60, 50);
end
