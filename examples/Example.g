
  #load the necessary packages
LoadPackage("GBNP");
LoadPackage("QPA");

  #define a quiver  Q with 4 vertices and 2 arrows from each vertex to the next
Q := Quiver(4, [ [1,2, "a1"], [1,2,"a2"], [2,3,"b1"], [2,3,"b2"], [3,4,"c1"], [3,4,"c2"], [4,1,"d1"], [4,1, "d2"]]);

  #This quiver comes from the ribbon graph ...

  #define the path algebra of Q over the field of rationals
kQ := PathAlgebra(Rationals, Q);

  #get a list of the generators of the path algebra.  one element for each vertex and one for each arrow of Q.
gens := GeneratorsOfAlgebra(kQ);

  #Rename the arrows (for convenience)
a1 := gens[5]; a2 := gens[6]; b1 := gens[7]; b2 := gens[8]; c1 := gens[9]; c2 := gens[10]; d1 := gens[11]; d2 := gens[12];

  #Make a list of relations
  #We start with one relation for each arrow.  These come from the super-potential, i.e., the faces
a1rel := b1*c1*d1 - b2*c1*d2;
a2rel := b2*c2*d2 - b1*c2*d1;
b1rel := c1*d1*a1 - c2*d1*a2;
b2rel := c2*d2*a2 - c1*d2*a1;
c1rel := d1*a1*b1 - d2*a1*b2;
c2rel := d2*a2*b2 - d1*a2*b1;
d1rel := a1*b1*c1 - a2*b1*c2;
d2rel := a2*b2*c2 - a1*b2*c1;

  #now add zigzag relations.  These are minimal paths (all of length 3, here) that are not part of any face.
zigzags := [a1*b1*c2, a1*b2*c2, a2*b2*c1, a2*b1*c1, b1*c1*d2, b1*c2*d2, b2*c2*d1, b2*c1*d1, c1*d1*a2, c1*d2*a2,
            c2*d2*a1, c2*d1*a1, d1*a1*b2, d1*a2*b2, d2*a2*b1, d2*a1*b1];

relations := [a1rel,a2rel,b1rel,b2rel,c1rel,c2rel,d1rel,d2rel];
Append(relations,zigzags);

  #Define the ideal of relations and form quotient PathAlgebra
I := Ideal(kQ, relations);

gb := GBNPGroebnerBasis(relations, kQ);

I := Ideal(kQ,gb);
GroebnerBasis(I,gb);

B := kQ/I;

  #Some commands to try
Dimension(B);
IsSymmetricAlgebra(B);
CartanMatrix(B);
IndecProjectiveModules(B);
