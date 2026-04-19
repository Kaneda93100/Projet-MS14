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
  write_sol_img("Joc_brut", Joc_brut);
  
  //--- Dummy Compression
  int comp_dummy = 0;

  if(comp_dummy == 1)
  { 
    printf("\nDébut de la compression avec dummy\n");
    to = clock();
    img* Joc_comp = dummy_comp(Joc_brut, 87);
    ti = clock();
    write_sol_img("Joc_comp_dummy", Joc_comp);
    double bazar = psnr(Joc_brut, Joc_comp);
    printf("Pnsr = %f\n", bazar);
    printf("RCT : %f\n", (double)(Joc_comp->M->NbrTri)/Joc_brut->M->NbrTri);
    printf("RCN : %f\n", (double)(Joc_comp->M->NbrVer)/Joc_brut->M->NbrVer);
    printf("Compression effectuée en %f sec\n", (ti - to)/CLOCKS_PER_SEC);    
  }

  //--- Random Compression
  int comp_random = 0;

  if(comp_random == 1)
  {
    printf("\nDébut de la compression aléatoire\n");
    to = clock();
    img* Joc_comp = random_comp(Joc_brut,0.1);
    ti = clock();
    write_sol_img("Joc_comp_dummy", Joc_comp);
    double bazar = psnr(Joc_brut, Joc_comp);
    printf("Pnsr = %f\n", bazar);
    printf("RCT : %f\n", (double)(Joc_comp->M->NbrTri)/Joc_brut->M->NbrTri);
    printf("RCN : %f\n", (double)(Joc_comp->M->NbrVer)/Joc_brut->M->NbrVer);
    printf("Compression effectuée en %f sec\n", (ti - to)/CLOCKS_PER_SEC);
    }

  
  //--- PSNR compression
  double target = 38;
  
  int comp_pnsr = 1;

  if(comp_pnsr == 1)
  {
    printf("\nDébut de la compression avec PNSR\n");
    printf("PNSR Target = %f", target);
    to = clock();
    img* Joc_comp = psnr_comp(Joc_brut, target);
    ti = clock();
    write_sol_img("Joc_comp_dummy", Joc_comp);
    double bazar = psnr(Joc_brut, Joc_comp);
    printf("Pnsr = %f\n", bazar);
    printf("RCT : %f\n", (double)(Joc_comp->M->NbrTri)/Joc_brut->M->NbrTri);
    printf("RCN : %f\n", (double)(Joc_comp->M->NbrVer)/Joc_brut->M->NbrVer);
    printf("Compression effectuée en %f sec\n", (ti - to)/CLOCKS_PER_SEC);    
  }

  int comp_pnsr_opt = 0;
  if(comp_pnsr_opt == 1)
  {
    printf("\nDébut de la compression avec PNSR optimisé\n");
    printf("PNSR Target = %f", target);
    to = clock();
    img* Joc_comp = psnr_comp_opt(Joc_brut, target);
    ti = clock();
    write_sol_img("Joc_comp_dummy", Joc_comp);
    double bazar = psnr(Joc_brut, Joc_comp);
    printf("Pnsr = %f\n", bazar);
    printf("RCT : %f\n", (double)(Joc_comp->M->NbrTri)/Joc_brut->M->NbrTri);
    printf("RCN : %f\n", (double)(Joc_comp->M->NbrVer)/Joc_brut->M->NbrVer);
    printf("Compression effectuée en %f sec\n", (ti - to)/CLOCKS_PER_SEC);    
  }

  return 0;
}