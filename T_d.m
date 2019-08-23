function value = T_d( u )

%% U_D evaluates the Dirichlet boundary conditions.
%
%  Discussion:
%
%    The user must supply the appropriate routine for a given problem
%
%  Modified:
%
%    23 February 2004
%
%  Parameters:
%
%    Input, real U(N,M), contains the M-dimensional coordinates of N points.
%
%    Output, VALUE(N), contains the value of the Dirichlet boundary
%    condition at each point.
%

  value=zeros(size(u,1),1);
  %valores para Ux
  for i=1:size(value,1)
  if(abs(u(i,1)-0)<0.001)
  value(i) = 1;
  else
  value(i) = 0;
  end
  end


  
