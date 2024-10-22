function pres = FunConstrainedOptimization(L,param_fiber,Deflection,Intersection,n_nodes,nb_stacks,lambda1,lambda2,eps)

parfevalOnAll(@warning,0,'off','all') 

tic
maxi = size(Deflection,1); %Number of fibers segmented
for j = 1:maxi 

    parfor t = 1:nb_stacks
        if size(Deflection{j,t}) == [n_nodes,3]
            j
            t
            intersection = Intersection{j,t};
            n_dof = 6;
            nd = n_nodes*n_dof;
            p = rand(nd,1)*10^-11; %generates a random force vector for the first iteration
            
            %Sets to 0 x and theta dof
            for i = 0:n_nodes-1 %forces and moments on degrees of freedom other than y and z are set to 0
                    p(i*n_dof+1) = 0;
                    p(i*n_dof+4) = 0;
                    p(i*n_dof+5) = 0;
                    p(i*n_dof+6) = 0;
            end
            %Sets to 0 nodes with no cell intersection
            for i = find(~intersection)
                    p((i-1)*n_dof+2) = 0;
                    p((i-1)*n_dof+3) = 0;
            end
            
            f = @(p)Calculate_cost_with_regularization_Constrained(L,param_fiber,n_nodes,p,intersection,Deflection{j,t},lambda1,lambda2,eps);           
            options = optimoptions('fminunc','Algorithm','quasi-newton','Diagnostics','on','SpecifyObjectiveGradient',true,'Display','iter','StepTolerance',1e-14);
            [pmin,fval] = fminunc(f,p,options); %optimization function
            
            pres{j,t} = reshape(pmin,[n_dof,n_nodes]);
            pres{j,t} = pres{j,t}'; %pres{j,t} contains the inferred force vector for fiber idx j at timepoint t
        end 
    end
end
toc
end

