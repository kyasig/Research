#Quiver has 2 vertices, each with in-deg = out-deg = 3, 2 loops at each vertex
# Bipartite graph has 2 6-faces and 2 vertices.
#face algebra is infinite-dimensional?

Read("../lib/lib.g");

n := 6;
e := epsilon(2*n);
f := (1,10,11,8,9,12)(2,5,4,7,6,3);
v := e*f;

G := getGraph(n,v); #Ribbon Graph Object
A := faceAlgebra(G);

Print("The dimension of A is ", Dimension(A), "\n");

Print("The Cartan matrix of A is \n", CartanMatrix(A), "\n");
P:=IndecProjectiveModules(A);
#Print("The indecomposable projective A-modules are \n", P, "\n");
Print("The radical series of P1 is \n", RadicalSeries(P[1]), "\n");

  #checking if A is symmetric can take a long time.
Print("Is A symmetric? ", IsSymmetricAlgebra(A));
