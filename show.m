function show ( elements3, coordinates, u )

%% SHOW displays the solution of the finite element computation.
%
%  Modified:
%
%    23 February 2004
%
%  Parameters:
%
%    Input, integer ELEMENTS3(N3,3), the nodes that make up each triangle.
%
%    Input, integer ELEMENTS4(N4,4), the nodes that make up each quadrilateral.
%
%    Input, real COORDINATES(N,1:2), the coordinates of each node.
%
%    Input, real U(N), the finite element coefficients which represent the solution.
%    There is one coefficient associated with each node.
%

%figure; clf

hold on
  triplot3 ( elements3, coordinates(:,1), coordinates(:,2), u');
%hold off

title ( 'Solution to the Problem' )
