#include "mesh.h"

int main(int argc, char* argv[])
{
  double to, ti;
  
  char* joc_msh = "../data/joconde.lowres.mesh";
  char* joc_sol = "../data/joconde.lowres.sol";

  //--- read a mesh
  to        = clock();
  Mesh* Joc_msh = msh_read(joc_msh, 1);
  double* Joc_sol = sol_read(joc_sol, Joc_msh->Dim, Joc_msh->NbrVer);
  ti        = clock();

  printf("  Vertices   %10d \n", Joc_msh->NbrVer);
  printf("  Triangles  %10d \n", Joc_msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);
  printf("\n\n\n");

  img* Joc_brut = init_from_mesh(Joc_msh, Joc_sol);
  img* Joc_comp = comp_img(Joc_brut);
  if (!Joc_brut)
    return 0;

  int succeed = msh_neighbors(Joc_comp->M);
  write_sol_img("Joc_comp", Joc_comp);
  
  return 0;
}