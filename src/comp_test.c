#include "mesh.h"

int main(int argc, char* argv[])
{
  double to, ti;
  
  char* joc_msh = "../data/joconde.mesh";
  char* joc_sol = "../data/joconde.sol";

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
  printf("\nDébut de la compression\n");
  img* Joc_precomp1 = dummy_comp_img(Joc_brut, 100);
  //img* Joc_precomp2 = dummy_comp_img(Joc_precomp1, 10);

  int print_precomp = 1;
  if(print_precomp == 1)
  {
    msh_neighbors(Joc_precomp1->M);
    write_sol_img("Joc_comp", Joc_precomp1);
    return 0;
  }
  img* Joc_comp = psnr_comp(Joc_precomp1, 35);
  if (!Joc_brut)
    return 0;

  msh_neighbors(Joc_comp->M);
  write_sol_img("Joc_comp", Joc_comp);
  printf("Jusqu'ici tout vas bien");
  printf("\nPSNR : %f\n", psnr(Joc_brut, Joc_comp));


  return 0;
}