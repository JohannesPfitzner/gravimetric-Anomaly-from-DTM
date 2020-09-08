function Gz = fTopographicReductionMagranaso(xObs,yObs,zObs,triangles,points,density)
%fTopographicReduction
%                                              
    X       = xObs;
    Y       = yObs;
    Z       = zObs;             
    Normals = 0;
    
    Gc = 6.6732e-6;             % Universal Gravitational constant
    [Nf,~] = size(triangles);   % Number of triangless
    [Ncor,~] = size(points);    % Number of Vertices
    Nedges=Nf*3;                % Number of Edges
    
    Edge=zeros(Nedges,8);       % initialize array
    Un=zeros(Nf,3);             % initialize array
    k=1;
    
    if isscalar(Normals)==0     %#ok
       Un=Normals;
       for f=1:Nf
       p1=points(triangles(f,1),:);
       p2=points(triangles(f,2),:);
       p3=points(triangles(f,3),:);
       sn=cross(p2-p1,p3-p1);
       UnTmp(f,:)=sn./norm(sn); %#ok
       if UnTmp(f,:)~=Un(f,:)
          idn(k)=f; %#ok
          k=k+1;
       end
       end
       tmp=triangles(idn,3);
       triangles(idn,3)=triangles(idn,1);
       triangles(idn,1)=tmp;
    else
       for f=1:Nf
       p1=points(triangles(f,1),:);
       p2=points(triangles(f,2),:);
       p3=points(triangles(f,3),:);
       sn=cross(p2-p1,p3-p1);
       Un(f,:)=sn./norm(sn);
       end
    end

    for f=1:Nf                  % Get edgelengths
       indx=[triangles(f,:) triangles(f,1)];
       for t=1:3
          edgeno=(f-1)*3+t;
          ends=indx(t:t+1);
          p1=points(ends(1),:);
          p2=points(ends(2),:);
          V=p2-p1;
          L=norm(V);
          Edge(edgeno,1:3)=V;
          Edge(edgeno,4)=L;
          Edge(edgeno,7:8)=ends;
       end
    end   
    
    [npro, nstn]=size(X);        % get mesh size  
    
    Gx=zeros(size(X)); 
    Gy=Gx; 
    Gz=Gx;
    
    % Comments: Now, for each observation point do the following: 
    % For each triangles find solid angle; 
    % for each side find p,q,r, and add p,q,r of sides 
    % to get P,Q,R for the triangles; 
    % find gx,gy,gz. 
    % Add the components from all the triangless to get 
    % Gx,Gy,Gz at the station.
    tic
    h=waitbar(0,'Calculating ...');
    for pr=1:npro
       for st=1:nstn
          opt=[X(pr,st) Y(pr,st) Z(pr,st)]; % GEÄNDERT Z
          fsign=zeros(1,Nf); 
          Omega=zeros(1,Nf); 
          for t=1:Ncor
             cor(t,:) = points(t,:)-opt; %#ok
          end % shift origin 

          for f=1:Nf
             nsides=3;
             cors=triangles(f,:); 
             Edge(:,5:6)=zeros(Nedges,2); % Clear record of integration 
             indx=[1:nsides 1 2];
             for t=1:nsides
                crs(t,:)=cor(cors(t),:); %#ok
             end
             % Find if the triangles is seen from inside
             fsign(f)=sign((Un(f,1)*crs(1,1)+Un(f,2)*crs(1,2)+Un(f,3)*crs(1,3)));
             % Find solid angle W subtended by triangles f at opt
             dp1=crs(indx(1),1)*Un(f,1)+crs(indx(1),2)*Un(f,2)+crs(indx(1),3)*Un(f,3);
             dp=abs(dp1); 
             if dp==0 
                Omega(f)=0;
             end
             if dp~=0
                
                tmp1=TriAngle(crs(3,:),crs(2,:),crs(1,:)); 

                Omega(f)=tmp1;
                if fsign(f)==-1
                    Omega(f)=-fsign(f)*Omega(f);
                end

             end 
             %indx=[1:nsides 1 2];
             %for t=1:nsides
             %   crs(t,:)=cor(cors(t),:);
             %end
             % Integrate over each side, if not done, and save result 
             PQR=[0 0 0];
             for t=1:nsides
                p1=crs(indx(t),:); 
                p2=crs(indx(t+1),:); 
                Eno=(f-1)*3+t; % Edge number 
                if Edge(Eno,6)==1
                   I=Edge(Eno,5);
                   V=Edge(Eno,1:3);
                   pqr=I .* V;
                   PQR=PQR+pqr; 
                end 
                if Edge(Eno,6)~=1  % in original manuscript there was a typo...
                   chsgn=1; % if origin, p1 & p2 are on a st line
                   if (p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3))/(norm(p1)*norm(p2))==1
                      if norm(p1)>norm(p2) % and pi farther than p2 
                         chsgn=-1; 
                         psave=p1; 
                         p1=p2; 
                         p2=psave; %#ok interchange p1,p2
                      end
                   end
                   V=Edge(Eno,1:3);
                   L=Edge(Eno,4);
                   L2=L*L; 
                   b=2*(V(1)*p1(1)+V(2)*p1(2)+V(3)*p1(3));
                   r1=norm(p1);
                   r12=r1*r1; 
                   b2=b/L/2; 
                   if (r1+b2 == 0)
                      V= -Edge(Eno,1:3);
                      b=2*(V(1)*p1(1)+V(2)*p1(2)+V(3)*p1(3));
                      b2=b/L/2; 
                   end 
                   if (r1+b2 ~= 0)
                      I = (1/L).* log((sqrt(L2 + b + r12) + L + b2)./(r1 + b2));
                   end 
                   s=find((Edge(:,7)==Edge(Eno,8)) & (Edge(:,8) == Edge(Eno,7)));
                   I=I*chsgn; % change sign of I if p1,p2 were interchanged 
                   Edge(Eno,5)=I;
                   Edge(s,5)=I;
                   Edge(Eno,6)=1;
                   Edge(s,6)=1;
                   pqr=I .* V; 
                   PQR=PQR+pqr; 
                end 
             end % From Omega,l,m,n,PQR, get components of field due to 
             % triangles f
             l=Un(f,1);
             m=Un(f,2);
             n=Un(f,3);
             p=PQR(1,1);
             q=PQR(1,2);
             r=PQR(1,3);

            if dp~=0 % if distance to triangles is non-zero
               gx=-density*Gc*dp1*(l*Omega(f)+n*q-m*r);
               Gx(pr,st)=Gx(pr,st)+ gx;
               gy=-density*Gc*dp1*(m*Omega(f)+l*r-n* p);
               Gy(pr,st)=Gy(pr,st)+ gy;
               gz=-density*Gc*dp1*(n*Omega(f)+m*p-l*q);
               Gz(pr,st)=Gz(pr,st)+ gz;
            end
          end
       end
       waitbar(pr/npro,h);
    end % end of triangless, stns, profiles
    close(h);
    toc
    
    Gz = -Gz;
    
end