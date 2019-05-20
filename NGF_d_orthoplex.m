function [a,kn] = NGF_d_orthoplex(N,s,beta,d,figure_l)

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
% with d-dimensional ORTHOPLEXES

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
% energy of the nodes epsilon is uniform from 0-9% Code that generates NGF in dimension d and flavour s=-1,0,1.
% (c) Ginestra Bianconi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization

     epsilon=floor(10*rand(N+2*d,1));
    
for n1=1:d,
    for s1=0:1,
        i1=(n1-1)*2+s1+1;
        for n2=1:(d),
            for s2=0:1,
                if(abs(n2-n1)>0),
                i2=(n2-1)*2+s2+1;
                a(i1,i2)=1;
                a(i2,i1)=1;
                end
            end
         end
        end
    end


for nt=1:(2^d),
    at(nt)=1;
    a_occ(nt)=1;
    a_occ3(nt)=1;
    v=de2bi(nt-1,d);
    for i=1:d,
        j=(i-1)*2+v(i)+1;
        node(nt,i)=j;
        at(nt)=at(nt)*exp(-beta*epsilon(j));        
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At each time t=in-8 we attach a new d-cube




it=2*d;
while it <=(N-d),
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
    
    
    %Add all the faces

    for n1=1:d,
        I2((n1-1)*2+1)=node(nj,n1);
        I2((n1-1)*2+2)=it+n1;
    end
    it=it+d;
    
        a(I2,I2)=a([1:(2*d)],[1:(2*d)]);
        
    for nt2=2:2^(d),
        nt=nt+1;
     at(nt)=1;
    a_occ(nt)=1;
    a_occ3(nt)=1;
    v=de2bi(nt2-1,d);
    for i=1:d,
        j=I2((i-1)*2+v(i)+1);
        node(nt,i)=j;
        at(nt)=at(nt)*exp(-beta*epsilon(j));        
    end
end
    end
k=sum(a>0);
v=2*(d-1);
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
