Read("../lib/lib.g");
#all 5-edge examples that are reduced without 2-faces have 1 10-face
#the quiver is always one vertex with 5 loops
#the dimension is 32 for each one, probably with radical layer dimensions 1,5,10,10,5,1.
#all appear to by symmetric (only checked first 4)
#Question: are any/all of these algebras isomorphic?
#Dimension of Center: 12, 9, 9, 7, 9, 9, 7, 9 => at least 3 distinct algebras up to isomorphism.
#Betti numbers for F1 and F2 agree (initially): 1,5,30,171, ? (takes too long)

n := 5;
v0 := (1,3,5,7,9);
  # edge perm is always (1,2)(3,4)(5,6)(7,8)(9,10)
  # vertex perm for Fk is v0*vk

v1 := (2,4,6,8,10);
g1 := getGraph(n,v0*v1);
F1 := faceAlgebra(g1);

Print("F1 has dimension ", Dimension(F1), ", with center of dimension ", Dimension(Center(F1)), "\n");
#Print("Is F1 symmetric? ", IsSymmetricAlgebra(F1), "\n");

S := SimpleModules(F1)[1];



v2 := (2,4,8,10,6);
g2 := getGraph(n,v0*v2);
F2 := faceAlgebra(g2);

Print("F2 has dimension ", Dimension(F2), ", with center of dimension ", Dimension(Center(F2)), "\n");
#Print("Is F2 symmetric? ", IsSymmetricAlgebra(F2), "\n");




v3 := (2,4,10,6,8);
g3 := getGraph(n,v0*v3);
F3 := faceAlgebra(g3);

Print("F3 has dimension ", Dimension(F3), ", with center of dimension ", Dimension(Center(F3)), "\n");


v4 := (2,6,10,4,8);
g4 := getGraph(n,v0*v4);
F4 := faceAlgebra(g4);

Print("F4 has dimension ", Dimension(F4), ", with center of dimension ", Dimension(Center(F4)), "\n");
#Print("Is F4 symmetric? ", IsSymmetricAlgebra(F4), "\n");



v5 := (2,6,8,4,10);
g5 := getGraph(n,v0*v5);
F5 := faceAlgebra(g5);

Print("F5 has dimension ", Dimension(F5), ", with center of dimension ", Dimension(Center(F5)), "\n");


v6 := (2,8,4,6,10);
g6 := getGraph(n,v0*v6);
F6 := faceAlgebra(g6);

Print("F6 has dimension ", Dimension(F6), ", with center of dimension ", Dimension(Center(F6)), "\n");




v7 := (2,8,4,10,6);
g7 := getGraph(n,v0*v7);
F7 := faceAlgebra(g7);

Print("F7 has dimension ", Dimension(F7), ", with center of dimension ", Dimension(Center(F7)), "\n");



v8 := (2,8,10,4,6);
g8 := getGraph(n,v0*v8);
F8 := faceAlgebra(g8);

Print("F8 has dimension ", Dimension(F8), ", with center of dimension ", Dimension(Center(F8)), "\n");


#AlgebraList := [F1,F2,F3,F4,F5,F6,F7,F8];

#for F in AlgebraList
