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
  return [Q,edgesOfQuiver,indices];
end;

getAlgebra := function(q)
 return PathAlgebra(Rationals,q);
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
  return Set(edges);
end;

#There doesnt seem to be a way to get the number of vertices, hence using the orignal graph
  #There should be. try "NumberOfVertices(q)", but there should be something using kq.
genArrows := function(kq,g)
  local arrowGens, n;
  arrowGens := GeneratorsOfAlgebra(kq);
  n := Length(Orbits(Group(getFaces(g)),[1..numEdges(g)]));
  return arrowGens{[n+1..Length(arrowGens)]};
end;

#no built in way to get the identity
  #have you tried "One(kq)" ?
identity := function(kq,g)
  local gens,sum, n;
  gens := GeneratorsOfAlgebra(kq);
  n := Length(Orbits(Group(getFaces(g)),[1..numEdges(g)]));
  sum := Sum(gens{[1..n]});
  return sum;
end;

multiplyEdges := function(halfEdges,q,kq,g)
 local edge, i,prod, arrows,ind;
 prod := identity(kq,g);
 arrows := genArrows(kq,g);
 for i in halfEdges do
  if i mod 2 = 0 then
    ind := i - 1;
  else
    ind := i;
  fi;
  prod := prod * arrows[ Position(q[3],ind)];
 od;
 return prod;
end;

superpotentialRelations := function(paths,q,kq,g)
  local relations,path,i,j, prod1,prod2,curr,next;
  relations := [];
  for i in [1..Length(paths)] do;
    if i mod 2 = 0 then
      continue;
    fi;
    curr := List(paths[i], x -> x[1]);
    next := List(paths[i+1], x -> x[1]);
    prod1 := multiplyEdges(curr{[2..Length(curr)]},q,kq,g);
    prod2 := multiplyEdges(next{[2..Length(next)]},q,kq,g);
    Add(relations, prod1-prod2);
  od;
  return relations;
end;

zigzagRelations := function(g)
  local edgeList,edge,phi,paths,fst;
  paths := [];
  phi := getFaces(g);
  epsilon := getEdges(g);
  edgeList := Orbits(Group(getEdges(g)), [1..numEdges(g)]);
  for edge in edgeList do
    fst := edge[1];
    Add(paths,[fst, ((fst +1) ^(phi^-1)), (fst ^ (epsilon * phi * epsilon * phi))]);
    Add(paths,[fst, ((fst+1)^phi), (fst ^ (epsilon * phi * epsilon * (phi^-1)))]);
  od;
  return paths;
end;

totalRelations := function(g,q,kq)
  local relations, superpotentialPaths,zigzags, i,j,k,halfEdges;
  relations := [];
  halfEdges := [];
  zigzags := zigzagRelations(g);
  for i in zigzags do
    Add(relations,multiplyEdges(i,q,kq,g));
  od;
  superpotentialPaths := superpotentialRelations(getSuperpotentialPaths(g),q,kq,g);
  Append(relations,superpotentialPaths);
  return relations;
end;

faceAlgebra := function(g)
  local i,q,kq, B;
  q:= buildQuiver(g);
  kq := getAlgebra(q[1]);
  i := Ideal(kq,totalRelations(g,q,kq));
  return kq/i;
  # may be able to replace "kq/i" by "kq/totalRelations(g,q,kq)" without defining the ideal, and this may even do the groebner basis too (?)
end;

