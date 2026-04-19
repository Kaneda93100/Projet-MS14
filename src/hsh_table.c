#include "mesh.h"

int main(int argc, char* argv[])
{
  char* mesh = "../data/squarecircle.mesh";

  int    iTri, iVer;
  double to, ti;
  
  //--- read a mesh
  to        = clock();
  Mesh* Msh = msh_read(mesh, 1);
  ti        = clock();

  //--- Tester les performances de msh_neighborsQ2 contre msh_neighbors
  to = clock();
  msh_neighborsQ2(Msh);
  ti = clock();

  printf("  time q2 neigh.        %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  free(Msh);
  Msh = NULL;

  to        = clock();
  Msh = msh_read(mesh, 1);
  ti        = clock();

  to = clock();
  msh_neighbors(Msh);
  ti = clock();

  printf("  time hsh method neigh.        %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //--- Test de la fonction de coloriage du domaine

  coloriage_magique(Msh);
  for(int i = 1; i <= Msh->NbrTri; i++)
  {
    printf("\n%d\n", Msh->TriMrk[i]);
  }
  msh_write(Msh, "colored_mesh.mesh");
  printf("\nMaillage colorié.\n");
  
  double* Qal = (double*)malloc(sizeof(double) * (Msh->NbrTri + 1));
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    Qal[iTri] = qual1(Msh,iTri);
    //printf("%d : %f\n", iTri, Qal[iTri]);
  }
  msh_write2dfield_Triangles("colored_mesh.solb", Msh->NbrTri, Qal);

  double3d* Met = (double3d*)malloc(sizeof(double3d) * (Msh->NbrVer + 1));  
  for (iVer = 1; iVer <= Msh->NbrVer; iVer++) {
    Met[iVer][0] = 1.;
    Met[iVer][1] = 0.;
    Met[iVer][2] = 1.;
  }
  msh_write2dmetric("metric.solb", Msh->NbrVer, Met);
  //--- Diverses informations sur la table de hachage
  printf("\n\n\n");
  
  HashTable* hsh = Hash_build(Msh);
  double av = Av_colision(hsh);
  printf("Nombre moyen de colision dans la table : %f", av);
  
  int* data = malloc(sizeof(int)*2);
  data = hash_biggest_head(hsh);

  printf("\n\nClé ayant le plus de colision : %d\nIl s'agit de la tête %d\n", data[0], data[1]);

  printf("\n\n\n");

  int2d* succeed = Edges_build(Msh);
  

  printf("\n\n\n");

  return 0;
}
