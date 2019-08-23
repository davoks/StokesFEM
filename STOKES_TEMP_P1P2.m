%---------------------------------------------------
% resolucion de las ecuaciones de Stokes en un dominio 2D 
% para resolver este problema elegimos elementos P2 para velocidades
% y elementos P1 para presiones y temperatura

  clear all;

  % cargamos los datos de la malla que los usaremos para resolver las presiones
  load coordinates.dat;
  load elements3.dat;
  load dirichlet.dat;

  % definimos una nueva malla para implementar las funciones P2 para las velocidades
  [nodes2element,nodes2edge,noedges,edge2element,interioredge]=edge(elements3,coordinates);

  elements3_p2=zeros(size(elements3,1),6);
  dirichlet_p2=zeros(2*size(dirichlet,1),2); 

  NUM_NODOS_P1=size(coordinates,1);
  NUM_NODOS_P2=noedges+NUM_NODOS_P1;
  coordinates_p2=zeros(NUM_NODOS_P2,2);
  nodos=[1:NUM_NODOS_P2];

  for i=1:NUM_NODOS_P1
      coordinates_p2(i,:)=coordinates(i,:);
  end

  for e=1:size(elements3,1);
      elements3_p2(e,[1 2 3])=elements3(e,:);
      n1=NUM_NODOS_P1+nodes2edge(elements3_p2(e,1),elements3_p2(e,2));
      elements3_p2(e,4)=n1;
      if(nodos(n1)!=0)
         coordinates_p2(n1,:)=(coordinates(elements3_p2(e,1),:)+coordinates(elements3_p2(e,2),:))/2;
         nodos(n1)=0;
      end
      n2=NUM_NODOS_P1+nodes2edge(elements3_p2(e,2),elements3_p2(e,3));
      elements3_p2(e,5)=n2;
      if(nodos(n2)!=0)
         coordinates_p2(n2,:)=(coordinates(elements3_p2(e,2),:)+coordinates(elements3_p2(e,3),:))/2;
         nodos(n2)=0;
      end
      n3=NUM_NODOS_P1+nodes2edge(elements3_p2(e,3),elements3_p2(e,1));
      elements3_p2(e,6)=n3;
      if(nodos(n3)!=0)
         coordinates_p2(n3,:)=(coordinates(elements3_p2(e,3),:)+coordinates(elements3_p2(e,1),:))/2;
         nodos(n3)=0;
      end
  end

  for i=1:size(dirichlet,1)
      dirichlet_p2(2*i-1,1)=dirichlet(i,1);
      dirichlet_p2(2*i-1,2)=NUM_NODOS_P1+nodes2edge(dirichlet(i,1),dirichlet(i,2));
      dirichlet_p2(2*i,1)=NUM_NODOS_P1+nodes2edge(dirichlet(i,1),dirichlet(i,2));
      dirichlet_p2(2*i,2)=dirichlet(i,2);
  end
%--------------------------------------------------------------------------------------------------------------
  % condiciones de Dirichlet para temperatura (x=0 y x=1 únicamente + flujo nulo en y=0 y y=1 
  % que es lo mismo que no poner condición de Dirichlet)
  diricount=1;
  dirichlet_T=zeros(size(dirichlet,1),2);
  for i=1:size(dirichlet,1)
      if((abs(coordinates(dirichlet(i,1),1)-1)<0.001 && abs(coordinates(dirichlet(i,2),1)-1)<0.001) ||
         (abs(coordinates(dirichlet(i,1),1)-0)<0.001 && abs(coordinates(dirichlet(i,2),1)-0)<0.001))
        dirichlet_T(diricount,:)=dirichlet(i,:);
	diricount+=1;
      end
  end


%--------------------------------------------------------------------------------------------------------------
printf("Problema de Stokes de dimensión %d \n",size(coordinates,2));
printf("Número de nodos P1= %d \n",NUM_NODOS_P1);
printf("Número de nodos P2= %d \n",NUM_NODOS_P2);
printf("Número de elementos= %d \n",size(elements3,1));
NUM_ITERACIONES=3;
X_n=zeros(2*NUM_NODOS_P2 + 2*NUM_NODOS_P1,1);
%for i=1:size(X_n,1)
%    X_n(i)=i;
%end
        BoundNodes_P2 = unique ( dirichlet_p2 );
        BoundNodes_P1 = unique ( dirichlet_T([1:diricount-1],:) );
        FreeNodes = setdiff ( 1:2*NUM_NODOS_P2+2*NUM_NODOS_P1, [BoundNodes_P2; BoundNodes_P2+NUM_NODOS_P2; 2*NUM_NODOS_P2+1 ; BoundNodes_P1+2*NUM_NODOS_P2+NUM_NODOS_P1]' );
        X_n([BoundNodes_P2 , BoundNodes_P2+NUM_NODOS_P2])   = u_d( [coordinates_p2(BoundNodes_P2,:) ; coordinates_p2(BoundNodes_P2,:)]);
        X_n([BoundNodes_P1+2*NUM_NODOS_P2+NUM_NODOS_P1])    = T_d( [coordinates(BoundNodes_P1,:)]);
        X_n(2*NUM_NODOS_P2+1,1)=1;
% Aplicamos Newton-Raphson para resolver el problema
for i=1:NUM_ITERACIONES
        [JAC,RES]=assStokesGuido(coordinates_p2, elements3_p2, X_n, NUM_NODOS_P1);
        X_n(FreeNodes) =X_n(FreeNodes)- JAC(FreeNodes,FreeNodes)\RES(FreeNodes);
        norm(RES(FreeNodes))
end

Ux=zeros(NUM_NODOS_P2,1);
Ux=X_n([1:NUM_NODOS_P2]);
Uy=zeros(NUM_NODOS_P2,1);
Uy=X_n([NUM_NODOS_P2+1 : 2*NUM_NODOS_P2]);
P=zeros(NUM_NODOS_P1,1);
P=X_n([2*NUM_NODOS_P2+1 : 2*NUM_NODOS_P2+NUM_NODOS_P1]);
T=zeros(NUM_NODOS_P1,1);
T=X_n([2*NUM_NODOS_P2+NUM_NODOS_P1+1 : 2*NUM_NODOS_P2+2*NUM_NODOS_P1]);

%figure;
%quiver(coordinates(:,1),coordinates(:,2),Ux([1:NUM_NODOS_P1]),Uy([1:NUM_NODOS_P1]),1)
%figure;
%show ( elements3, coordinates, full ( P ) );
%figure;
%show ( elements3, coordinates, full ( T ) );

save -ascii "temp.txt" coordinates T

%----------------------------------------------------------------------------------------------  

