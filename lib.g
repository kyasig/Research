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

get_matrix_directed := function(n,p)
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
  local vertexList, edgeList, adjMatrix, s,t, n,i, Q;
  n := g[4];
  edgeList := [];
  vertexList := Orbits(Group(getFaces(g)),[1..n]);
  adjMatrix := get_matrix_directed(n, getFaces(g));
  for s in [1..Length(adjMatrix)] do
    for t in [1..Length(adjMatrix)] do
      if adjMatrix[s][t] = 0 then
        continue;
      fi;
      for i in [1..(adjMatrix[s][t]) ] do
        Add(edgeList, [s,t,"a"]);
      od;
    od;
  od;
  Q := Quiver(Length(vertexList), edgeList);
  return Q;
end;
