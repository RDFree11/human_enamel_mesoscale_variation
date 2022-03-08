function [XYZv,conn,XYp] = intersect_plane_surf(q,nh,Xg,Yg,Zg,DT,tol)
% intersect_plane_surf determines the intersection of a plane and a
% triangulated surface
% the plane is define by point q and unit normal vector nh
% the surface is defined by meshgrids Zg, Yg, and Zg = F(Xg,Yg)
% DT is the 2D Delaunay triangulation of (Xg,Yg)
% tol is a threshhold to identify close-by points in the intersection

    % calculate signed distance of all vertices from plane
    XYZv = [Xg(:),Yg(:),Zg(:)];
    d = (XYZv-q')*nh;
    
    % create array of signed distances for vertices of all triangles
    s = sign(d(DT));
    
    % determine the number of triangle vertices that fall onto the plane
    % todo: implemement tolerance criterion
    n = zeros(length(s),1);
    for ii = 1:length(s)
        n(ii) = 3-nnz(s(ii,:));
    end
    
    % initialize connectivity matrix
    % each row describes contains two indices of vertices in XYZv
    % the intersection is a line between the indexed vertices
    conn = [];
    
    % initialize counter for vertices added later
    ctr = length(XYZv)+1;
    
    % if all three vertices are on the plane, add 3 edges to connectivity
    % matrix
    ind3 = find(n==3);
    for ii = 1:length(ind3)
        V = DT(ind3(ii),s(ind3(ii),:)==0);
        conn = [conn;V(1:2);V([1,3]);V(2:3)];
    end
    
    % if exactly two vertices are on the plane, add 1 edge to connectivity matrix
    ind2 = find(n==2);
    for ii = 1:length(ind2)
        conn = [conn;DT(ind2(ii),s(ind2(ii),:)==0)];
    end
    
    % if one vertex is on the plane
    ind1 = find(n==1);
    for ii = 1:length(ind1)
        % and if the other two vertices are on opposite sides, find intercept of plane
        % with crossing edge
        if sum(nonzeros(s(ind1(ii),:)))==0
            % two vertices, q and r, that are on opposite side of the plane
            qr = XYZv(DT(ind1(ii),s(ind1(ii),:)~=0),:);
            % the vertex p that is on the plane
            p  = XYZv(DT(ind1(ii),s(ind1(ii),:)==0),:);
            % unit normal to the plane through the vertices
            nth = cross(p-qr(1,:),p-qr(2,:));
            nth = nth/norm(nth);
            % the unit direction vector of the line of intersection between
            % the section plane and the triangle plane
            dh1 = cross(nth,nh);
            % the unit direction vector of the triangle edge that crosses
            % the plane
            dh2 = diff(qr)/norm(diff(qr));
            % solve for intercept between the line through p with direction
            % vector dh1 and the the qr edge
            dv = [dh1;dh2]'\(qr(1,:)-p)';
            % add intercept to vertex array and increment counter
            XYZv(ctr,:) = dv(1)*dh1+p;
            ctr = ctr+1;
            % update connectivity matrix
            conn = [conn;[DT(ind1(ii),s(ind1(ii),:)==0),ctr]];           
        end
        % if the other vertices are on the same side, ignore
    end
    
    % find all remaining triangles (w/o vertices on the
    % plane) that have vertices on opposite sides of the plane
    % note that there are the following possibilities for a row in s:
    % +/-[1,1,1]; all permutations of +/-[-1,1,1]
    % if we take the absolute of the row sum, we either get 3 or 1
    % only the latter have vertices on opposite sides
    indt = find(n == 0 & abs(sum(s,2))==1);
    
    % iterate over all triangles with edges that cross the plane
    for ii = 1:length(indt)
        % grab vertices
        V = XYZv(DT(indt(ii),:),:);
        % get row vector from s for triangle
        sv = s(indt(ii),:);
        % make sure that we have a permutation of [1,-1,-1] 
        if sum(sv)==1
            sv = -1*sv;
        end
        % now vertex p is on ones side
        p  = V(sv==1,:);
        % vertices q and r on the other
        qr = V(sv==-1,:);
        
        % for the two edges that cross the plane, find the intercept of the
        % edge with the plane
        for jj=1:2
            % unit vector of direction of crossing edge
            dh = (-p+qr(jj,:))/norm(-p+qr(jj,:));
            % solve for intersection between edge and plane
            d  = ((q'-p)*nh)/(dh*nh);
            % determine new vertex of line segment
            XYZv(ctr,:) = d * dh + p;
            ctr = ctr+1;
        end
        % update connectivity matrix
        conn = [conn;[ctr-2,ctr-1]];
    end
    
    % identify unique vertices and their number
    uverts = unique(conn);
    Nv     = length(uverts);
    % determine distance between vertices
    D = nan(Nv);
    for ii=1:Nv
        for jj=ii+1:Nv
            D(jj,ii) = norm(XYZv(uverts(ii),:)-XYZv(uverts(jj),:));
        end
    end
    % find vertices that are closer than a threshhold value
    [jind,iind] = find(D <= tol);
    
    
    for ii = 1:length(jind)
        % replace one of the close-by vertices with their mean
        XYZv(uverts(iind(ii)),:) = mean(XYZv([uverts(iind(ii)),uverts(jind(ii))],:));
        % and replace reference to the other vertex in the connectivity list 
        conn(conn == uverts(jind(ii))) = uverts(iind(ii));
    end
    % to identify connected regions, determine the graph
    G = graph(categorical(conn(:,1)),categorical(conn(:,2)));
    % identify connected components of the graph
    [bins,binsize] = conncomp(G);
    Ng = unique(bins);
    
    for ii = 1:length(Ng)
        sG = subgraph(G,bins == ii);
        beg = find(degree(sG)==1,1,'first');
        
        if isempty(beg)
            beg = 1;
        end
        sorted_node_indices = dfsearch(sG,beg);
        sorted_indices = str2num(char(sG.Nodes{sorted_node_indices,1}));
            
        % convert to basis of plane   
        if (~nh(1)) & ~(nh(2)) & (nh(3)) % if the plane normal is parallel to the z-axis
            XYZp = [[0;0;sign(nh(3))],[1;0;0],[0;1;0]]^-1*(XYZv(sorted_indices,:)'-q);
        else 
            b1 = nh;
            b2 = cross([0;0;1],b1);
            b2 = b2/norm(b2);
            b3 = cross(b1,b2);
            XYZp = [b1,b2,b3]^-1*(XYZv(sorted_indices,:)'-q);
        end
        XYp{ii} = XYZp(2:3,:);
    end
end

