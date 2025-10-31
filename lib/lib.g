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

get_matrix_directed := function(n,p) #TODO make this just take a graph g as inpute instead
  local a,e,pair,i,j,source,target,v, ai,aj;
  e := epsilon(n);
  v := Orbits(Group(p),[1..n]);
  a := NullMat(Length(v),Length(v));
  for pair in Orbits(Group(e),[1..n]) do
    i := pair[1];
    j := pair[2];
    for source in v do
      if i in source then
        ai := Position(v,source);
        for target in v do
         if j in target then
           aj := Position(v,target);
           a[ai][aj] := a[ai][aj] + 1;
         fi;
        od;
      fi;
    od;
  od;
  return a;
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
           name := [CharInt(96 + source), CharInt(48+ind)]; #trick for converting to chars
            #it might not be necessary to name the arrows -- just omit the "name" argument. But it doesn't hurt.
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


#There doesnt seem to be a way to get the number of vertices, hence using the orignal graph
  #There should be. try "NumberOfVertices(q)", but there should be something using kq.
genArrows := function(kq)
  local arrowGens, n;
  arrowGens := GeneratorsOfAlgebra(kq);
  n := NumberOfVertices(QuiverOfPathAlgebra(kq));
  return arrowGens{[n+1..Length(arrowGens)]};
end;


multiplyEdges := function(halfEdges,kq,indices)
 local edge, i,prod, arrows,ind;
 prod := One(kq);
 arrows := genArrows(kq);
 for i in halfEdges do
  if i mod 2 = 0 then
    ind := i - 1;
  else
    ind := i;
  fi;
  prod := prod * arrows[Position(indices,ind)];
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
  local hedgeList,edge,phi,paths,fst;
  paths := [];
  phi := getFaces(g);
  epsilon := getEdges(g);
  hedgeList := Orbits(Group(getEdges(g)), [1..numEdges(g)]);
  for edge in hedgeList do
    fst := edge[1];
    Add(paths,[fst, ((fst +1) ^(phi^-1)), (fst ^ (epsilon * phi * epsilon * phi))]);

    Add(paths,[fst, ((fst+1)^phi), (fst ^ (epsilon * phi * epsilon * (phi^-1)))]);
  od;
  return paths;
end;

zigzagRelations := function(paths,kq,indices)
  local relations, i;
  relations := [];
  for i in paths do
    Add(relations,multiplyEdges(i,kq,indices));
  od;
  return relations;
end;

totalRelations := function(g,kq,indices)
  local relations,zigzags;
  relations := superpotentialRelations(getSuperpotentialPaths(g),kq,indices);
  zigzags := zigzagRelations(zigzagPaths(g),kq,indices);
  Append(relations,zigzags);
  return Set(relations);
end;

faceAlgebra := function(g)
  local i,q,kq;
  q := buildQuiver(g);
  kq := PathAlgebra(Rationals,q[1]);
  return kq/totalRelations(g,kq,q[2]);
  # may be able to replace "kq/i" by "kq/totalRelations(g,q,kq)" without defining the ideal, and this may even do the groebner basis too (?)
end;

