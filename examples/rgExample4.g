#Quiver has 5 vertices, each with in-deg = out-deg = 2.
# Bipartite graph has 5 4-faces and 3 vertices.

Read("../lib/lib.g");

n := 10;
e := epsilon(2*n);
f := (1,18,11,10)(2,3,20,13)(4,5,12,15)(6,7,14,17)(8,9,16,19);
v := e*f;

G := getGraph(10,v); #Ribbon Graph Object
A := faceAlgebra(G);

Print("The dimension of A is ", Dimension(A), "\n");

Print("The Cartan matrix of A is \n", CartanMatrix(A), "\n");
P:=IndecProjectiveModules(A);
#Print("The indecomposable projective A-modules are \n", P, "\n");
Print("The radical series of P1 is \n", RadicalSeries(P[1]), "\n");

  #checking if A is symmetric can take a long time.
#Print("Is A symmetric? ", IsSymmetricAlgebra(A));
