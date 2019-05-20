function [a,kn] = NGF_d_hypercube(N,s,beta,d,figure_l)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  If you use this code, please cite:
%  G. Bianconi and C. Rahmede 
%  "Network geometry with flavour: from complexity to quantum geometry"
%  Physical Review E 93, 032315 (2016)
%  and
%  D. Mulder and G. Bianconi
% "Network geometry and complexity"
%  Journal of Statistical Physics 173,783 (2018).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code that generates NGF in dimension d and flavour s=-1,0,1 costructed
% with d-dimensional HYPERCUBES

% a adjacency matrix
% kn vector of  generalized degrees   of the nodes


% This code uses 
% N maximal number of nodes in the NGF
% Flavour of the NGF  s=-1,0,1
% Inverse temperature: beta>0 or beta=0
% Dimension d with d>1
% figure_l=1 will print the edge list of the network in file 
%"NGF_edgelist_d%d_s%d.edges"
% figure_l=0 will not print the edge list of the network
% energy of the nodes epsilon is uniform from 0-9
% (c) Ginestra Bianconi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign energies to the nodes
% If using Poisson and power-law you must define
% the parameters mu, or kappa
% Examples:
% mu=10;
% kappa=1;

     epsilon=floor(10*rand(N+2^d,1));
% Alternative energy distributions
    %epsilon(i)=random('Poisson',mu); 
    %poisson distribution with average mu
    %epsilon(i)=rand(1)^(1/(kappa+1)); 
    %power-law distribution with exponent kappa
   



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Initial condition: at time t=1 a single d-dimensional hypercube 
for j1=1:(2^d),
    v1=de2bi(j1-1,d);
    for j2=1:(2^d),       
        v2=de2bi(j2-1,d);
        if nnz((v1-v2))==1,
            a(j1,j2)=1;
            a(j2,j1)=1;
        end
    end
end
        
  nt=0;

  
for m1=1:d,
    for np1=0:1,
        nt=nt+1;
        n1=0;
        at(nt)=1;
        a_occ(nt)=1;
        a_occ3(nt)=1;
        for j=1:(2^d),
           v=de2bi(j-1,d);
           if (v(m1)==np1),
               n1=n1+1;
               node(nt,n1)=j;
               at(nt)=at(nt)*exp(-beta*epsilon(j));
           end
       end
    end
end
            
        
    
it=2^d;

while it<(N-2^(d-1))
    
% Choose (d-1) face to which to attach the new d-cube

    [I,J,V]=find(at.*a_occ);

    norm=sum(V);
    x=rand(1)*norm;
    for nj1=1:numel(V),
            x=x-V(nj1);
         if x<0,
             nj=J(nj1);
             break;
         end
    end    
    
     a_occ(nj)=a_occ(nj)+s; %update of 1+sn
     a_occ3(nj)=a_occ3(nj)+1;%update of the generalized degree
     
     
     %In=zeros(2^d,1);
     for n2=1:2^(d-1)
         In(n2)=node(nj,n2);
         In(2^(d-1)+n2)=it+n2;
     end
    
     a(In,In)=a([1:2^(d)],[1:(2^d)]);
    it=it+2^(d-1);
     for m1=1:d,    
    for np1=0:1,
        if(((d-m1+np1)>0))
            nt=nt+1;
        n1=0;
        at(nt)=1;
        a_occ(nt)=1;
        a_occ3(nt)=1;
        for j=1:2^d,
           v=de2bi(j-1,d);
           if((v(m1)==np1)),
               n1=n1+1;
               node(nt,n1)=In(j);
               at(nt)=at(nt)*exp(-beta*epsilon(In(j)));
           end
        end
       end
        end
    end
end
        
     
     
     
    k=sum(a>0);
    k=sum(a>0);
    v=d;
    kn=(k-v)./(v-d+1)+1;

 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% % Print network file
 [I2,J2,A2]=find(tril(a>0));
 if (figure_l==1)
    filename=sprintf('NGF_edgelist_d%d_s%d.edges',d,s);
    fid=fopen(filename,'w');
    for it=1:numel(A2),
        fprintf(fid, ' %d %d  \n', I2(it), J2(it));
    end
    fclose(fid);
 end
end
