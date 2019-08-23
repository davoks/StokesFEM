%
%  Esta función ensambla la matriz Jacobiana y el Residuo en una malla de simplices n-D
%  para resolver el problema de Stokes con elementos P1 Y P2 para (vx,vy) Y ,p0 y T
%
%  Argumentos:
%  * coordinates: matriz de coordenadas de la malla
%  (nnodos x dimensiones)
%  * elements: nodos de cada elemento
%  (nelementos x (dimensiones + 1))
%  * degree: grado del polinomio a integrar exactamente, determina
%  la cuadratura a utilizar.
%  * X_n solución en el paso anterior
%
function [JAC, RES] = assStokesGuido (coordinates_p2, elements3_p2, X_n, NUM_NODOS_P1)

g=9.8;
betta=100;
mu=1;
k=0.01;

NUM_NODOS_P2=size(coordinates_p2,1);
weights=load("2D/strang4_w.txt");
gpoints=load("2D/strang4_x.txt");
% Definimos JAC y RES
JAC=sparse(2*NUM_NODOS_P2+2*NUM_NODOS_P1,2*NUM_NODOS_P2+2*NUM_NODOS_P1);
RES=sparse(2*NUM_NODOS_P2+2*NUM_NODOS_P1,1);

% Definimos los vectores que tendrán el valor de las funciones de forma evaluados en los puntos de Gauss
phi_p2   = zeros(6,1); % 6 funciones de forma por triángulo
d_phi_p2 = zeros(6,2);          
phi_p1   = zeros(3,1); % 3 funciones de forma por triángulo
d_phi_p1 = zeros(3,2);             

% Recorremos cada elemento
for e = 1 : size(elements3_p2,1) 

        % matrices elementales y vectores elementales
	f1_v = zeros(6,6); 
        f1_p = zeros(6,3);
	f2_p = zeros(6,3);
        f1_t = zeros(6,3); 
        f4_vx= zeros(3,6);
        f4_vy= zeros(3,6);
        f4_t = zeros(3,3);
	res_1= zeros(6,1); 
        res_2= zeros(6,1);
        res_3= zeros(3,1);
        res_4= zeros(3,1);

        for xg = 1 : size(gpoints,1)

            for i=1:3
        	phi_p1(i)      =   PHI_2D(i,1,gpoints(xg,:)); 
        	d_phi_p1(i,:)  = D_PHI_2D(i,1,gpoints(xg,:));
	    endfor
            for i=1:6
       		phi_p2(i)      =   PHI_2D(i,2,gpoints(xg,:)); 
		d_phi_p2(i,:)  = D_PHI_2D(i,2,gpoints(xg,:));
	    endfor
            % calculamos el Jacobiano de la transformación lineal para calcular el valor del gradiente 
            % correctamente en el elemento en que estamos
            JACOBI=coordinates_p2(elements3_p2(e,[1 2 3]),:)' * d_phi_p1;
 	    %transformamos las coordenadas del punto de gauss en el master al punto en el elemento real
            x_e   =coordinates_p2(elements3_p2(e,[1 2 3]),:)' * phi_p1;
            M=det(JACOBI)*weights(xg);

                % calculamos las variables en el punto de gauss
                % estás serviran para construir los términos no lineales
                % en el JAC y para el RES
                ux=X_n(elements3_p2(e,:))'*phi_p2;
                uy=X_n(elements3_p2(e,:)+NUM_NODOS_P2)'*phi_p2;
                P =X_n(elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2)'*phi_p1;
                T =X_n(elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2+NUM_NODOS_P1)'*phi_p1;
                dux=[0;0];
                duy=[0;0];
                dP =[0;0];
                dT =[0;0];
                for i=1:6
                	dux+=X_n(elements3_p2(e,i))*(JACOBI'\d_phi_p2(i,:)');
                	duy+=X_n(elements3_p2(e,i)+NUM_NODOS_P2)*(JACOBI'\d_phi_p2(i,:)');
                end
                for i=1:3
                	dP +=X_n(elements3_p2(e,i)+2*NUM_NODOS_P2)*(JACOBI'\d_phi_p1(i,:)');
                	dT +=X_n(elements3_p2(e,i)+2*NUM_NODOS_P2+NUM_NODOS_P1)*(JACOBI'\d_phi_p1(i,:)');
                end
                %---------------------------------------------------
        	% construimos las matrices elementales
                for  i = 1 : 6
           		for  j = 1 : 6
                		f1_v(i,j) += (JACOBI'\d_phi_p2(i,:)')' * (JACOBI'\d_phi_p2(j,:)')*M;         
           		end      
                        res_1(i,1)+= ((JACOBI'\d_phi_p2(i,:)')'*dux-P*(JACOBI'\d_phi_p2(i,:)')'*[1;0])*M;
                        res_2(i,1)+= ((JACOBI'\d_phi_p2(i,:)')'*duy-P*(JACOBI'\d_phi_p2(i,:)')'*[0;1]-g*betta*T*phi_p2(i))*M;            
      		end

        	for  i = 1 : 6
           		for  j = 1 : 3
                		f1_p(i,j) += (JACOBI'\d_phi_p2(i,:)')'*[1;0] * phi_p1(j)*M ;
				f2_p(i,j) += (JACOBI'\d_phi_p2(i,:)')'*[0;1] * phi_p1(j)*M ;
                                f1_t(i,j) += phi_p2(i)*phi_p1(j)                        *M ;
                                f4_vx(j,i)+= phi_p1(j)*phi_p2(i)*dT(1)                  *M ;
                                f4_vy(j,i)+= phi_p1(j)*phi_p2(i)*dT(2)                  *M ;
           		endfor                   
      		endfor

                for  i = 1 : 3
                     for  j = 1 : 3
                         f4_t(i,j) +=([ux uy]*(JACOBI'\d_phi_p1(j,:)')*phi_p1(i) +k*(JACOBI'\d_phi_p1(i,:)')' * (JACOBI'\d_phi_p1(j,:)')) *M;
                     end
                     res_3(i,1)+= -( dux(1)+duy(2) )*phi_p1(i)                            *M ;
                     res_4(i,1)+= ([ux uy]*dT*phi_p1(i) + k*(JACOBI'\d_phi_p1(i,:)')'*dT) *M ;
                end

        end % end Gauss points

        %Ensamblamos JAC y RES
        %-------------------
        % vx | 0  |  p  | 0 
        % 0  | vy |  p  | T
        % vx | vy |  0  | 0
        % vx | vy |  0  | T 
        %-------------------        


        JAC(elements3_p2(e,:) , elements3_p2(e,:))                                                                      += mu*f1_v;
        JAC(elements3_p2(e,:)+NUM_NODOS_P2, elements3_p2(e,:)+NUM_NODOS_P2)                                             += mu*f1_v;
        JAC([elements3_p2(e,:)],[elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2])                                               += -f1_p;
        JAC([elements3_p2(e,:)+NUM_NODOS_P2],[elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2])                                  += -f2_p;
        JAC([elements3_p2(e,:)+NUM_NODOS_P2],[elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2+NUM_NODOS_P1])                     += -g*betta*f1_t;
        %-------------------
        % vx | 0  |  p  | 0 
        % 0  | vy |  p  | T
        %-------------------        

        JAC(elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2,elements3_p2(e,:))                                                   += -f1_p';
        JAC(elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2,elements3_p2(e,:)+NUM_NODOS_P2)                                      += -f2_p';
        % vx | vy |  0  | 0 

        JAC(elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2+NUM_NODOS_P1,elements3_p2(e,:))                                      +=f4_vx;
        JAC(elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2+NUM_NODOS_P1,elements3_p2(e,:)+NUM_NODOS_P2)                         +=f4_vy;
        JAC(elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2+NUM_NODOS_P1,elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2+NUM_NODOS_P1)    +=f4_t;
        % vx | vy |  0  | T 
        %-------------------         

        RES(elements3_p2(e,:),1)                                                                                        +=res_1;
        RES(elements3_p2(e,:)+NUM_NODOS_P2,1)                                                                           +=res_2;
        RES(elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2,1)                                                                   +=res_3;
        RES(elements3_p2(e,[1 2 3])+2*NUM_NODOS_P2+NUM_NODOS_P1,1)                                                      +=res_4;


end

endfunction

