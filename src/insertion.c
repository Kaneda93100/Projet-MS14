#include "mesh.h"

int main()
{
  /*
    Dans ce script, on teste la fonction insertion autour de laquelle
    la majorité du code est construite. 
    La variable vers_ins est le nombre de point que l'on souhaite insérer
    dans le maillage "qube", qui est un carré de côté 1. 
  */

  int ver_ins = 1000;
  Mesh* qube = init_dummy_qube(ver_ins);

  if (!qube)
  {
    printf("\nMesh not found :(\n");
    return 0;
  }
  printf("\nMesh founded !\n");

  double G[2];
  G[0] = 0; G[1] = 0;


  srand(42); // La réponse à la grande question

  double start_insert = clock();
  for(int i = 0; i < ver_ins; i++)
  {
    // Générer un point dans le carré
    G[0] = ((double)rand())/RAND_MAX;
    G[1] = ((double)rand())/RAND_MAX;

    // Insérer le point G dans le qube
    insertion(qube, G);
    printf("\n%d/%d\n\n%d-eme ajout\n\n", qube->NbrVer, qube->NbrVerMax, i+1);
  }

  msh_write(qube, "qube_debile.mesh");
  double stop_insert = clock();

  return 0;
}