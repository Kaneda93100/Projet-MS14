#include "mesh.h"

int main(int argc, char* argv[])
{

  printf("Bad mozerfucker\n");

  
  int    iTri, iVer;
  double to, ti;
  
  //--- read a mesh
  to        = clock();
  Mesh* Msh = msh_read(argv[1], 1);
  ti        = clock();

  if (!Msh)
    return 0;

  printf("  Vertices   %10d \n", Msh->NbrVer);
  printf("  Triangles  %10d \n", Msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  printf("\n\n\n");

  //--- create neigbhors with hash table
  int succeed = msh_neighbors(Msh);
  volatile HashTable* hsh = Hash_build(Msh);

  for(int i = 0; i < hsh->SizHead; i++){printf("%d\t", hsh->Head[i]);}

  double av = Av_colision(hsh);
  printf("Nombre moyen de colision dans la table : %f", av);
  
  int* data = malloc(sizeof(int)*2);
  data = hash_biggest_head(hsh);

  printf("\n\nClé ayant le plus de colision : %d\nIl s'agit de la tête %d\n", data[0], data[1]);

  printf("\n\n\n");

  succeed = Edges_build(Msh);

  Edges_vertices_cout(Msh);

  printf("\n\n\n");

  //--- create neigbhors Q2 version
  // to = clock();
  // msh_neighborsQ2(Msh);
  // ti = clock();

  // printf("  time q2 neigh.        %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //--- TODO: compute mesh quality
  double* Qal = (double*)malloc(sizeof(double) * (Msh->NbrTri + 1));



  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    Qal[iTri] = (double)iTri / 10.;
  }

  msh_write2dfield_Triangles("quality.solb", Msh->NbrTri, Qal);

  //--- TODO: compute metric field
  double3d* Met = (double3d*)malloc(sizeof(double3d) * (Msh->NbrVer + 1));

  //double Qual1 = qual1(Msh, 1);
  double Qual2 = qual2(Msh, 1);

  return 0;

  for (iVer = 1; iVer <= Msh->NbrVer; iVer++) {
    Met[iVer][0] = 1.;
    Met[iVer][1] = 0.;
    Met[iVer][2] = 1.;
  }

  msh_write2dmetric("metric.solb", Msh->NbrVer, Met);

  //--- Free memory
  if (Qal != NULL) {
    free(Qal);
    Qal = NULL;
  }
  if (Met != NULL) {
    free(Met);
    Met = NULL;
  }

  return 0;
}
