LoadPackage("QPA");
LoadPackage("GBNP");

epsilon := function(n)
  local e, i;
  e := ();
  for i in [1,3..n-1] do
    e := e * (i,i+1);
  od;
  return e;
end;

getGraph := function(n,p)
 local g, edges;
 edges := epsilon(2*n);
 g := [edges, p, (edges * p) , (2*n)];
 return g;
end;

getEdges := function(g)
  return g[1];
end;

getVertices := function(g)
  return g[2];
end;

getFaces := function(g)
  return g[3];
end;

numEdges := function(g)
  return g[4];
end;

buildQuiver := function(g)
  local vertexList, edgeList,edgesOfQuiver, Q,source,target,i, n, ind, name,indices;
  n := numEdges(g);
  ind := 1; #This is just for naming each arrow in the quiver
  vertexList := Orbits(Group(getFaces(g)),[1..n]);
  edgesOfQuiver := [];
  indices := []; #List of odd half-edges, their position in this list corresponds to the position of their associated edge in the quiver
  for source in [1..Length(vertexList)] do
    for i in vertexList[source] do
      if (i mod 2 = 1) then
        for target in [1..Length(vertexList)] do
          if i+1 in vertexList[target] then
           name := [CharInt(96 + source), CharInt(48+ind)];
           Add(edgesOfQuiver,[source,target, name]);
           ind := ind +1;
           Add(indices,i);
          fi;
        od;
      fi;
    od;
    ind := 1;
  od;
  Q := Quiver(Length(vertexList), edgesOfQuiver);
  return [Q,indices];
end;

getRightFace := function(edge, g)
  local curr, face;
  face := [edge[1]];
  curr := (edge[1]) ^ getVertices(g);
  while curr <> edge[1] do
    Add(face, curr);
    curr := curr ^ (getVertices(g));
  od;
  return face;
end;

getLeftFace := function(edge, g)
  local curr, face;
  face := [edge[2]];
  curr := (edge[2]) ^ (getVertices(g)^-1);
  while curr <> edge[2] do
    Add(face, curr);
    curr := curr ^ (getVertices(g)^-1);
  od;
  return face;
end;


genArrows := function(kq)
  local arrowGens, n;
  arrowGens := GeneratorsOfAlgebra(kq);
  n := NumberOfVertices(QuiverOfPathAlgebra(kq));
  return arrowGens{[n+1..Length(arrowGens)]};
end;

multiplyEdges := function(halfEdges,kq,indices) #uses a list of indices that match a half-edge with its position in the quiver
 local edge, i,prod, arrows,halfEdge;
 prod := One(kq);
 arrows := genArrows(kq);
 for halfEdge in halfEdges do
  if halfEdge mod 2 = 0 then
    i := halfEdge - 1;
  else
    i := halfEdge;
  fi;
  prod := prod * arrows[Position(indices,i)];
 od;
 return prod;
end;

getSuperpotentialPaths:= function(g)
  local edgeList, edges,right,left, i, j,index,pathList;
  edges := [];
  pathList := [];
  edgeList := Orbits(Group(getEdges(g)), [1..numEdges(g)]);
  for i in edgeList do
    right := getRightFace(i,g);
    for j in right do
      index := Int((j+1)/2);
      #Ceil throws some strange error for whatever reason
      Add(pathList, edgeList[index]);
    od;
    Add(edges, pathList);
    pathList := [];
    left := getLeftFace(i,g);
    for j in left do
      index := Int((j+1)/2);
      Add(pathList, edgeList[index]);
    od;
    Add(edges, pathList);
    pathList := [];
  od;
  return edges;
end;

superpotentialRelations := function(paths,kq,indices)
  local relations,path,i,j, prod1,prod2,curr,next;
  relations := [];
  for i in [1..Length(paths)] do;
    if i mod 2 = 0 then
      continue;
    fi;
    curr := List(paths[i], x -> x[1]);
    next := List(paths[i+1], x -> x[1]);
    prod1 := multiplyEdges(curr{[2..Length(curr)]},kq,indices);
    prod2 := multiplyEdges(next{[2..Length(next)]},kq,indices);
    Add(relations, prod1-prod2);
  od;
  return relations;
end;

zigzagPaths := function(g)
  local edgeList,edge,phi,paths,fst,snd,phiInv,epsilon;
  paths := [];
  phi := getFaces(g);
  phiInv := phi ^-1;
  epsilon := getEdges(g);
  edgeList := Orbits(Group(getEdges(g)), [1..numEdges(g)]);
  for edge in edgeList do
    fst := edge[1]; #2i-1
    snd := edge[2]; #2i
    Add(paths,
      [fst, snd ^ phiInv, (snd ^ ( phiInv * epsilon * phi))]
    );
    Add(paths,
      [fst, snd ^ phi, (snd ^ (phi * epsilon * phiInv))]
    );
  od;
  return paths;
end;

otherZigzagPaths:= function(g,q,indices)
  local edgeList,i,j,k,arrow,vertex,correspondingHedge,paths,fst,filtered,nu;
  edgeList := Orbits(Group(getEdges(g)), [1..numEdges(g)]);
  nu := getVertices(g);
  paths := [];
  filtered := [];
  for i in edgeList do
    fst := i[1]; #2i-1
    arrow := ArrowsOfQuiver(q)[Position(indices,fst)];
    vertex := TargetOfPath(arrow);
    if OutDegreeOfVertex(vertex) > 2 then
      for j in OutgoingArrowsOfVertex(vertex) do
        correspondingHedge := indices[Position(ArrowsOfQuiver(q),j)];
        Add(paths,[fst,correspondingHedge]);
      od;
    fi;
  od;
  for j in paths do
    if ((j[1]^nu) = j[2]) then #check if edge is part of a face
      Add(filtered,j);
    fi;
    if ((j[1]+1) = (j[2]+1)^nu) then #check if edge is part of a face
      Add(filtered,j);
    fi;
  od;
  return Difference(paths,filtered); #filter out edges part of a face
end;

zigzagRelations := function(paths,kq,indices)
  local relations, i;
  relations := [];
  for i in paths do
    Add(relations,multiplyEdges(i,kq,indices));
  od;
  return relations;
end;


totalRelations := function(g,kq,indices,q)
  local relations,zigzags,otherZigzags;
  relations := superpotentialRelations(getSuperpotentialPaths(g),kq,indices);
  zigzags := zigzagRelations(zigzagPaths(g),kq,indices);
  otherZigzags := zigzagRelations(otherZigzagPaths(g,q,indices),kq,indices);
  Append(relations,zigzags);
  Append(relations,otherZigzags);
  #relations := Filtered(relations,x -> not IsEmpty(x));
  return DuplicateFreeList(relations);
end;

faceAlgebra := function(g)
  local i,q,kq;
  q := buildQuiver(g);
  kq := PathAlgebra(Rationals,q[1]);
  return kq/totalRelations(g,kq,q[2],q[1]);
end;
