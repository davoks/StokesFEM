% está función calcula el valor de la función de forma phi_k de orden n
% en el elemento master en un punto xg en 2D
% está función fue creada para calcular los valores en los puntos
% de Gauss y así poder integrarlas numéricamente

function value = PHI_2D ( k , n , xg)

phi_1= det([1,xg(1,:);1 1 0;1 0 1])/det([1 0 0;1 1 0;1 0 1]);
phi_2= det([1,xg(1,:);1 0 0;1 0 1])/det([1 1 0;1 0 0;1 0 1]);
phi_3= det([1,xg(1,:);1 0 0;1 1 0])/det([1 0 1;1 0 0;1 1 0]);

x=xg(1);
y=xg(2);

if (n==1)
	if (k==1)
                value = phi_1;
		return; 
	endif
	if (k==2)
		value = phi_2;
		return;
	endif
	if (k==3)
		value = phi_3;
		return ;
	endif
elseif (n==2)
        if (k==1)
		value = (1-2*x-2*y)*(1-x-y);
		return;
	endif
	if (k==2)
		value = x*(2*x-1);
		return; 
	endif
	if (k==3)
		value = y*(2*y-1);
		return; 
	endif
	if (k==4)
		value = 4*x*(1-x-y);
		return; 
	endif
	if (k==5)
		value = 4*x*y;
		return;
	endif
	if (k==6)
	        value = 4*y*(1-x-y);
		return;
	endif
elseif (n==3)
        if (k==1)
		value = phi_1*(3*phi_1-1)*(3*phi_1-2)/2;
		return;
	endif
	if (k==2)
		value = phi_2*(3*phi_2-1)*(3*phi_2-2)/2;
		return; 
	endif
	if (k==3)
		value = phi_3*(3*phi_3-1)*(3*phi_3-2)/3;
		return; 
	endif
	if (k==4)
		value = 9*phi_1*phi_2*(3*phi_1-1)/4;
		return; 
	endif
	if (k==5)
		value = 9*phi_1*phi_2*(3*phi_2-1)/4;
		return;
	endif
	if (k==6)
	        value = 9*phi_3*phi_2*(3*phi_2-1)/4;
		return;
	endif
	if (k==7)
		value = 9*phi_3*phi_2*(3*phi_1-1)/4;
		return; 
	endif
	if (k==8)
		value = 9*phi_3*phi_1*(3*phi_3-1)/4;
		return;
	endif
	if (k==9)
	        value = 9*phi_3*phi_1*(3*phi_1-1)/4;
		return;
	endif
	if (k==10)
	        value = 27*phi_1*phi_2*phi_3;
		return;
	endif
endif
return;
end
