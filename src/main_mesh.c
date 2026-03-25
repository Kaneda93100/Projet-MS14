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

  // printf("\nAffichage de Msh->Efr\n");
  // for(int i = 1; i <= Msh->NbrEfrMax; i++){printf("%d : (%d, %d)\n",i,Msh->Efr[i][0], Msh->Efr[i][1]);}
  

  // printf("\nAffichage de Msh->EfrRef\n");
  // for(int i = 1; i <= Msh->NbrEfrMax; i++){printf("%d : %d\n", i, Msh->EfrRef[i]);}
  
  printf("  Vertices   %10d \n", Msh->NbrVer);
  printf("  Triangles  %10d \n", Msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  printf("\n\n\n");


  int succeed = msh_neighbors(Msh);
  printf("%d", succeed);

  // Tests des fonctions de calcul du rayon du cercle circonscrit et de son centre.
  double ray = circ_circle_ray(Msh, 1);
  printf("Rayon du cercle circonscrit du triangle [P0 = (%f,%f), P1 = (%f,%f), P2 = (%f,%f)]\nray = %f",Msh->Crd[Msh->Tri[1][0]][0], Msh->Crd[Msh->Tri[1][0]][1]
                                                                                                       ,Msh->Crd[Msh->Tri[1][1]][0], Msh->Crd[Msh->Tri[1][1]][1] 
                                                                                                       ,Msh->Crd[Msh->Tri[1][2]][0], Msh->Crd[Msh->Tri[1][2]][1]
                                                                                                       ,ray);

  double* centre = centre_circ_circle(Msh, 1);
  printf("\nCentre du cercle circonscrit de tri_equi : (%f,%f)\n", centre[0], centre[1]);
  
  if(succeed != 1){printf("Erreur dans le calcul des voisins.\n"); return 0;}


  double* P = malloc(sizeof(double)*2);
  P[0] = 0.45; P[1] = 0.125;
  
  // Test de la fonction de calcul de la cavité
  int2d* cavity = compute_cavity(Msh, 5, P);
  int edg = 0;
  while(cavity[edg][0] == 0 && cavity[edg][1] == 0){edg++;}

  for(int i = 0; i < edg; i++){printf("\n%d : (%d, %d)\n", i, cavity[i][0], cavity[i][1]);}

  for(int i = 1; i <= Msh->NbrTri; i++)
  {
    printf("\nTriangle %d : ", i);
    if(empty_sphere_criterion(Msh, i, P) == 1)
    {
      printf("1\n");
    }
    else
    {
      printf("0\n");
    }
  }
  
  return 0;
  // Test de la fonction d'insertion d'un point dans un maillage
  printf("\nMsh->NbrVerMax current : %d", Msh->NbrVerMax);
  push_NbrVerMax(Msh);
  printf("\nMsh->NbrVerMax updated : %d", Msh->NbrVerMax);


  insertion(Msh,P);
  for(int i = 1; i <= Msh->NbrVer; i++)
  {
    printf("\nNoeud %d : (%f,%f)\n", i, Msh->Crd[i][0], Msh->Crd[i][1]);
  }
  for(int i = 1; i <= Msh->NbrTriMax; i++)
  {
    printf("\nTriangle numéro %d : (%d,%d,%d)\n", i, Msh->Tri[i][0]
                                                   , Msh->Tri[i][1]
                                                   , Msh->Tri[i][2]);
  }

  // int delaunay_test_4 = empty_sphere_criterion(Msh, 4, P);
  // printf("Résultat pour itri == 3 : %d", delaunay_test_4);
  // double* alpha = malloc(sizeof(double)*2); alpha[0] = 0.65; alpha[1] = 0.15;
  // int pos = localising(Msh, alpha);
  // printf("Le point (%f,%f) est localisé dans le triangle %d\n", alpha[0], alpha[1], pos);
  // int delaunay_test_2 = empty_sphere_criterion(Msh, 2, alpha);
  // printf("Résultat pour itri == 2 : %d \n", delaunay_test_2);
  // double* bary_itri1 = malloc(sizeof(double)*3);
  // bary_itri1 = compute_BC(Msh, tri_ref, P);
  // for(int i = 0; i < 3; i++){printf("barry %d == %f\n", i, bary_itri1[i]);}
  // double * P1;
  // P1 = malloc(sizeof(double)*2); P1[0] = 0.25; P1[1] = 0.12;
  // int true = is_point_localised(Msh, tri_ref, P1);
  // int false = is_point_localised(Msh, tri_ref, P);
  // printf("\ntrue == %d\nfalse == %d\n", true, false);
  // int vert1 = 1, vert2 = 2;
  // int voi = identify_voi(Msh, tri_ref, vert1, vert2);
  // printf("vert1 : %d,\nvert2 : %d\n",vert1, vert2);
  // printf("\nvoi : %d", voi);
  // int loc_P = localising(Msh, P);
  // printf("Point P = (%f,%f) localisé dans le triangle %d\n", P[0], P[1],loc_P);

  return 0;


  coloriage_magique(Msh);
  for(int i = 1; i <= Msh->NbrTri; i++){printf("%d : %d\n",i,Msh->TriRef[i]);}
  msh_write(Msh, "colored.mesh");

  double* Qal = (double*)malloc(sizeof(double) * (Msh->NbrTri + 1));
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    Qal[iTri] = qual1(Msh,iTri);
    //printf("%d : %f\n", iTri, Qal[iTri]);
  }
  msh_write2dfield_Triangles("quality.solb", Msh->NbrTri, Qal);


  double3d* Met = (double3d*)malloc(sizeof(double3d) * (Msh->NbrVer + 1));
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
  // //--- create neigbhors with hash table
  // int succeed = msh_neighbors(Msh);
  // volatile HashTable* hsh = Hash_build(Msh);

  // // for(int i = 0; i < hsh->SizHead; i++){printf("%d\t", hsh->Head[i]);}

  // double av = Av_colision(hsh);
  // printf("Nombre moyen de colision dans la table : %f", av);
  
  // int* data = malloc(sizeof(int)*2);
  // data = hash_biggest_head(hsh);

  // printf("\n\nClé ayant le plus de colision : %d\nIl s'agit de la tête %d\n", data[0], data[1]);

  // printf("\n\n\n");

  // succeed = Edges_build(Msh);
  

  // printf("\n\n\n");

  //--- create neigbhors Q2 version
  to = clock();
  msh_neighborsQ2(Msh);
  ti = clock();

  printf("  time q2 neigh.        %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  return 0;
}
