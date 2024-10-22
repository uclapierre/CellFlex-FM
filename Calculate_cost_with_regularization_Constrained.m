function [Cost,Gradient]=Calculate_cost_with_regularization_Constrained(L,param_fiber,nb_nodes,p,intersection,deflection,lambda1,lambda2,eps)
   
pois = 0.5; rho_st = 2000.;
E_st = 11.32*1e6; G_st = E_st / (2*(1+pois));

 if param_fiber == 600
     x = 72.7*1e-9;
     z = 1.5567*1e-6;
     T = 7.05697517606636e-08;

 elseif param_fiber == 650
     x = 84*1e-9;
     z = 1.7567*1e-6;
     T = 1.27030199738801e-07;
 elseif param_fiber == 700
     x = 108.3*1e-9;
     z = 1.8933*1e-6;
     T = 1.79868456539690e-07;
 elseif param_fiber == 750
     x = 117.3*1e-9;
     z = 2.1800*1e-6;
     T = 2.69659548154013e-07;
 elseif param_fiber == 800
     x = 142*1e-9;
     z = 2.2333*1e-6;
     T = 3.73822606980116e-07;
 end

 A = x*z;
 Iy = (x*z^3)/12;
 Iz = (z*x^3)/12;
 J = Iy + Iz;
 Ep_poutre = [E_st G_st A Iy Iz J A*rho_st];

 % Model Topology
 ne = nb_nodes-1; % Number of elements
 Coord (1,:) = [0 0 0];
 for i = 1:ne
     Coord (1+i,:) = [ Coord(i,1)+L/ne ,0, 0 ];
     Elem(i,:) = [ i, i+1 ];
 end
 ePoutre = [1:ne];
 n_fix = [1,ne+1]; % Fixed nodes
 [n_nodes,n_dof,n_elem,n_nel,Dof,Edof] = topol (Coord,Elem);
 [Ex,Ey,Ez] = coordxtr(Edof,Coord,Dof,n_nel);

 %Matrices K et M
 nd = n_nodes*n_dof;
 K = zeros(nd); M = zeros(nd); C = zeros(nd);
 for ie = 1:ne
     eo(ie,:) = [0 0 1];
     [ke,me] = beam3eUnderTension(Ex(ie,:),Ey(ie,:),Ez(ie,:),eo(ie,:),Ep_poutre);
     K = assem(Edof(ie,:),K,ke);
     M = assem(Edof(ie,:),M,me);
 end
    
 %Boundary conditions
 bc = []; [b,bc,nb] = fix_point (bc, n_fix, Dof);
 %Forces and Moments calculation
 [~, ~, xyzF] = fe_stat (K,p,b,n_dof,n_nodes);

 %Optimization step - Cost calculation
 Diff = xyzF(:,2:3) - deflection(:,2:3);
 Cost = sum(abs(Diff(:))) + lambda1*(sum(abs(p(:)))) + lambda2*(sum(abs(p(:))))^2 ; %Cost with regularization (L1 : lambda1 + L2 : lambda2)

 %Optimization step - Gradient calculation
 dCost = zeros(size(p));
     
    for i = find(intersection)-1 %We calculate the gradient only for nodes with cell mask intersection, other nodes have gradient = 0

        k = i*n_dof+2; %i-th degree of freedom related to y deflection
        p_k = p;
        p_k(k) = p_k(k) + eps;
        [~, ~, xyzFk] = fe_stat (K,p_k,b,n_dof,n_nodes);
        Diff_k = xyzFk(:,2:3) - deflection(:,2:3);
        Cost_k = sum(abs(Diff_k(:))) + lambda2*(sum(abs(p_k(:))))^2+lambda1*(sum(abs(p_k(:)))) ;
        dCost(k) = Cost_k - Cost;
        
        k = i*n_dof+3; %i-th degree of freedom related to z deflection
        p_k = p;
        p_k(k) = p_k(k) + eps;
        [~, ~, xyzFk] = fe_stat (K,p_k,b,n_dof,n_nodes);
        Diff_k = xyzFk(:,2:3) - deflection(:,2:3);
        Cost_k = sum(abs(Diff_k(:))) + lambda2*(sum(abs(p_k(:))))^2+lambda1*(sum(abs(p_k(:)))) ;
        dCost(k) = Cost_k - Cost;
    end
    Gradient = dCost/eps; 
    
end