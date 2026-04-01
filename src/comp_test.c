#include "mesh.h"

int main(int argc, char* argv[])
{
  double to, ti;
  
  char * joc = "../data/joconde.lowres.mesh";

  //--- read a mesh
  to        = clock();
  Mesh* Joc_brut = msh_read(joc, 1);
  ti        = clock();

  if (!Joc_brut)
    return 0;
  
  printf("  Vertices   %10d \n", Joc_brut->NbrVer);
  printf("  Triangles  %10d \n", Joc_brut->NbrTri);
  printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  printf("\n\n\n");
 
  Mesh* Joc_comp = comp_img(Joc_brut, 5, bernoulli_criterion_comp);

  msh_write(Joc_comp, "test.mesh");

  return 0;
}