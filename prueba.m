% 
  phi_p2   = zeros(6,1); % 6 funciones de forma por triángulo
  d_phi_p2 = zeros(6,2);          
  phi_p1   = zeros(3,1); % 3 funciones de forma por triángulo
  d_phi_p1 = zeros(3,2);       
  gp=[0.0 0.25];
           for i=1:3
        	phi_p1(i)      =   PHI_2D(i,1,gp); 
        	d_phi_p1(i,:)  = D_PHI_2D(i,1,gp);
	    endfor
            for i=1:6
       		phi_p2(i)      =   PHI_2D(i,2,gp); 
		d_phi_p2(i,:)  = D_PHI_2D(i,2,gp);
	    endfor
            % calculamos el Jacobiano de la transformación lineal para calcular el valor del gradiente correctamente
            JACOBI=coordinates_p2(elements3_p2(3,[1 2 3]),:)' * d_phi_p1;
            d_phi_p1(1,:)'
            JACOBI\d_phi_p1(3,:)'
            JACOBI'\d_phi_p1(3,:)'

            JACOBI'\d_phi_p2(1,:)'

            det(JACOBI)

