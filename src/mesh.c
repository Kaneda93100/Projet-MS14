#include "mesh.h"

int tri2edg[3][2] = { { 1, 2 }, { 2, 0 }, { 0, 1 } };

Mesh* msh_init()
{
  Mesh* Msh = malloc(sizeof(Mesh));
  if (!Msh) return NULL;

  Msh->Dim    = 0;
  Msh->NbrVer = 0;
  Msh->NbrTri = 0;
  Msh->NbrEfr = 0;
  Msh->NbrEdg = 0;

  Msh->NbrVerMax = 0;
  Msh->NbrTriMax = 0;
  Msh->NbrEfrMax = 0;
  Msh->NbrEdgMax = 0;

  Msh->Box[0] = 1.e30; // xmin
  Msh->Box[1] = -1.e30; // xmax
  Msh->Box[2] = 1.e30; // ymin
  Msh->Box[3] = -1.e30; // ymax

  //--- Data for the list of vertices
  Msh->Crd = NULL;

  //--- Data for the list of triangles
  Msh->Tri    = NULL;
  Msh->TriVoi = NULL;
  Msh->TriRef = NULL;
  Msh->TriMrk = NULL;

  //--- Data for the list of boundary edges
  Msh->Efr    = NULL;
  Msh->EfrVoi = NULL;
  Msh->EfrRef = NULL;

  //--- Data for the list of edges
  Msh->Edg = NULL;

  return Msh;
}
Mesh* msh_read(char* file, int readEfr)
{
  char   InpFil[1024];
  float  bufFlt[2];
  double bufDbl[2];
  int    i, bufTri[4], bufEfr[3];
  int    FilVer, ref;

  int fmsh = 0;

  if (!file) return NULL;

  Mesh* Msh = msh_init();

  //--- set file name
  strcpy(InpFil, file);
  if (strstr(InpFil, ".mesh")) {
    if (!(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim))) {
      return NULL;
    }
  }
  else {
    strcat(InpFil, ".meshb");
    if (!(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim))) {
      strcpy(InpFil, file);
      strcat(InpFil, ".mesh");
      if (!(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim))) {
        return NULL;
      }
    }
  }

  printf(" File %s opened Dimension %d Version %d \n", InpFil, Msh->Dim, FilVer);

  Msh->NbrVer = GmfStatKwd(fmsh, GmfVertices);
  Msh->NbrTri = GmfStatKwd(fmsh, GmfTriangles);

  Msh->NbrVerMax = Msh->NbrVer;
  Msh->NbrTriMax = Msh->NbrTri;

  //--- allocate arrays
  Msh->Crd    = calloc((Msh->NbrVerMax + 1), sizeof(double3d));
  Msh->Tri    = calloc((Msh->NbrTriMax + 1), sizeof(int3d));
  Msh->TriRef = calloc((Msh->NbrTriMax + 1), sizeof(int1d));
  Msh->TriMrk = calloc((Msh->NbrTriMax + 1), sizeof(int1d));

  //--- read vertices
  GmfGotoKwd(fmsh, GmfVertices);
  if (Msh->Dim == 2) {
    if (FilVer == GmfFloat) { // read 32 bits float
      for (i = 1; i <= Msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufFlt[0], &bufFlt[1], &ref);
        Msh->Crd[i][0] = (double)bufFlt[0];
        Msh->Crd[i][1] = (double)bufFlt[1];
      }
    }
    else { // read 64 bits float
      for (i = 1; i <= Msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufDbl[0], &bufDbl[1], &ref);
        Msh->Crd[i][0] = bufDbl[0];
        Msh->Crd[i][1] = bufDbl[1];
      }
    }
  }
  else {
    fprintf(stderr, "  ## ERROR: 3D is not implemented\n");
    exit(1);
  }

  //--- read triangles
  GmfGotoKwd(fmsh, GmfTriangles);
  for (i = 1; i <= Msh->NbrTri; ++i) {
    GmfGetLin(fmsh, GmfTriangles, &bufTri[0], &bufTri[1], &bufTri[2], &bufTri[3]);
    Msh->Tri[i][0] = bufTri[0];
    Msh->Tri[i][1] = bufTri[1];
    Msh->Tri[i][2] = bufTri[2];
    Msh->TriRef[i] = bufTri[3];
  }

  //--- read boundary edges
  if (readEfr == 1) {
    Msh->NbrEfr    = GmfStatKwd(fmsh, GmfEdges);
    Msh->NbrEfrMax = Msh->NbrEfr;

    Msh->Efr    = calloc((Msh->NbrEfrMax + 1), sizeof(int2d));
    Msh->EfrRef = calloc((Msh->NbrEfrMax + 1), sizeof(int1d));

    GmfGotoKwd(fmsh, GmfEdges);
    for (i = 1; i <= Msh->NbrEfr; ++i) {
      GmfGetLin(fmsh, GmfEdges, &bufEfr[0], &bufEfr[1], &bufEfr[2]);
      Msh->Efr[i][0] = bufEfr[0];
      Msh->Efr[i][1] = bufEfr[1];
      Msh->EfrRef[i] = bufEfr[2];
    }
  }

  GmfCloseMesh(fmsh);

  return Msh;
}
double* sol_read(char* file, int mshDim, int mshNbrSol)
{
  char   InpFil[1024];
  int    FilVer, SolTyp, NbrTyp, SolSiz, TypTab[GmfMaxTyp];
  float  bufFlt;
  double bufDbl;
  int    i, dim, nbrSol;

  int fsol = 0;

  if (!file) return NULL;

  double* sol = NULL;

  //--- set file name
  strcpy(InpFil, file);
  if (strstr(InpFil, ".sol")) {
    if (!(fsol = GmfOpenMesh(InpFil, GmfRead, &FilVer, &dim))) {
      return NULL;
    }
  }
  else {
    strcat(InpFil, ".solb");
    if (!(fsol = GmfOpenMesh(InpFil, GmfRead, &FilVer, &dim))) {
      strcpy(InpFil, file);
      strcat(InpFil, ".sol");
      if (!(fsol = GmfOpenMesh(InpFil, GmfRead, &FilVer, &dim))) {
        return NULL;
      }
    }
  }

  printf(" File %s opened Dimension %d Version %d \n", InpFil, dim, FilVer);

  SolTyp = GmfSolAtVertices; // read only sol at vertices
  nbrSol = GmfStatKwd(fsol, SolTyp, &NbrTyp, &SolSiz, TypTab);

  if (nbrSol == 0) {
    printf("  ## WARNING: No SolAtVertices in the solution file !\n");
    return NULL;
  }
  if (dim != mshDim) {
    printf("  ## WARNING: WRONG DIMENSION NUMBER. IGNORED\n");
    return NULL;
  }
  if (nbrSol != mshNbrSol) {
    printf("  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    return NULL;
  }
  if (NbrTyp != 1) {
    printf("  ## WARNING: WRONG FIELD NUMBER. IGNORED\n");
    return NULL;
  }
  if (TypTab[0] != GmfSca) {
    printf("  ## WARNING: WRONG FIELD TYPE. IGNORED\n");
    return NULL;
  }

  sol = (double*)calloc(nbrSol + 1, sizeof(double));

  GmfGotoKwd(fsol, SolTyp);

  for (i = 1; i <= nbrSol; ++i) {
    if (FilVer == GmfFloat) {
      GmfGetLin(fsol, SolTyp, &bufFlt);
      sol[i] = (double)bufFlt;
    }
    else {
      GmfGetLin(fsol, SolTyp, &bufDbl);
      sol[i] = bufDbl;
    }
  }

  if (!GmfCloseMesh(fsol)) {
    fprintf(stderr, "  ## ERROR: Cannot close solution file %s ! \n", InpFil);
    // myexit(1);
  }

  return sol;
}
int msh_boundingbox(Mesh* Msh)
{
  int1d iVer;
  double Mx,My = 0;
  double mx, my = 0;

  //--- compute bounding box
  for (iVer = 1; iVer <= Msh->NbrVer; iVer++) {
      double tempx = Msh->Crd[iVer][0];
      double tempy = Msh->Crd[iVer][1]; 

      //Affectation du max
      if(Mx < tempx){Mx = tempx;}
      if(My < tempy){My = tempy;}

      //Affectation du min
      if(mx > tempx){mx = tempx;}
      if(my > tempy){my = tempy;}
  }

  Msh->Box[0] = Mx; 
  Msh->Box[1] = mx; 
  Msh->Box[2] = My; 
  Msh->Box[3] = my;


  return 1;
}
int msh_write(Mesh* Msh, char* file)
{
  int iVer, iTri, iEfr;
  int FilVer = 2;

  if (!Msh) return 0;
  if (!file) return 0;

  int fmsh = GmfOpenMesh(file, GmfWrite, FilVer, Msh->Dim);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  GmfSetKwd(fmsh, GmfVertices, Msh->NbrVer);
  for (iVer = 1; iVer <= Msh->NbrVer; iVer++)
    GmfSetLin(fmsh, GmfVertices, Msh->Crd[iVer][0], Msh->Crd[iVer][1], 0);

  GmfSetKwd(fmsh, GmfTriangles, Msh->NbrTri);
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++)
    GmfSetLin(fmsh, GmfTriangles, Msh->Tri[iTri][0], Msh->Tri[iTri][1], Msh->Tri[iTri][2], Msh->TriRef[iTri]);

  if (Msh->NbrEfr > 0) {
    GmfSetKwd(fmsh, GmfEdges, Msh->NbrEfr);
    for (iEfr = 1; iEfr <= Msh->NbrEfr; iEfr++)
      GmfSetLin(fmsh, GmfEdges, Msh->Efr[iEfr][0], Msh->Efr[iEfr][1], Msh->EfrRef[iEfr]);
  }

  GmfCloseMesh(fmsh);

  return 1;
}
int msh_neighborsQ2(Mesh* Msh)
{
  int iTri, iEdg, jTri, jEdg, iVer1, iVer2, jVer1, jVer2;

  if (!Msh) return 0;

  if (Msh->TriVoi == NULL)
    Msh->TriVoi = calloc((Msh->NbrTri + 1), sizeof(int3d));

  //--- Compute the neighbors using a quadratic-complexity algorithm
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    for (iEdg = 0; iEdg < 3; iEdg++) {
      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];

      //--- find the Tri different from iTri that has iVer1, iVer2 as vertices
      for (jTri = 1; jTri <= Msh->NbrTri; jTri++) {
        if (iTri == jTri)
          continue;

        for (jEdg = 0; jEdg < 3; jEdg++) {
          jVer1 = Msh->Tri[jTri][tri2edg[jEdg][0]];
          jVer2 = Msh->Tri[jTri][tri2edg[jEdg][1]];
          
          if((jVer1 == iVer1 && jVer2 == iVer2) || (jVer1 == iVer2 && jVer2 == iVer1))
          {
            Msh->TriVoi[iTri][iEdg] = jTri;
            Msh->TriVoi[jTri][jEdg] = iTri;
          }
          else
          {
            continue;
          }
        }
      }
    }
  }

  return 1;
}

// Méthodes sur les HashTable
volatile HashTable* hash_init(int SizHead, int NbrMaxObj) 
{
  HashTable* hsh = malloc(sizeof(HashTable)); 

  hsh->SizHead = SizHead; 
  hsh->NbrMaxObj = NbrMaxObj;
  hsh->NbrObj = 0;

  hsh->Head = malloc(sizeof(int)*SizHead);
  for(int i = 0; i < SizHead; i++){hsh->Head[i] = 0;}
  hsh->LstObj = malloc(sizeof(int5d)*NbrMaxObj);
  hsh->LstObj[0][0] = 0; hsh->LstObj[0][1] = 0; hsh->LstObj[0][2] = 0;hsh->LstObj[0][3] = 0; hsh->LstObj[0][4] = 0;

  return hsh;
}
int hash_find(volatile HashTable* hsh, int iVer1, int iVer2) 
{ 
  int k = hsh->Head[iVer1+iVer2];

  if(hsh->NbrObj == 0){
    //printf("Table vide.\n");
    return 0; 
  }

  while(k != 0){
    if((hsh->LstObj[k][0] == iVer1 && hsh->LstObj[k][1] == iVer2) || (hsh->LstObj[k][0] == iVer2 && hsh->LstObj[k][1] == iVer1)){return k;} // Arête trouvée
    else{k = hsh->LstObj[k][4];}
  }
  return 0; // Arête non trouvée
}
int hash_add(volatile HashTable* hsh, int iVer1, int iVer2, int iTri) 
{
  if(hsh->NbrObj + 1 > hsh->NbrMaxObj)
  {
    printf("La liste est complète. Pas d'ajout.\n");
    return 0;
  } 

  if(hsh->Head[iVer1+iVer2] !=0) // Tête non-vide (Les têtes sont toutes initialisées à 0)
  {
    int pos = hash_find(hsh, iVer1, iVer2); 
    if(pos == 0)
    {
      // Création du nouvel élément (placé en queu de liste)
      hsh->LstObj[hsh->NbrObj+1][0] = iVer1;
      hsh->LstObj[hsh->NbrObj+1][1] = iVer2;
      hsh->LstObj[hsh->NbrObj+1][2] = iTri;
      hsh->LstObj[hsh->NbrObj+1][3] = 0;
      hsh->LstObj[hsh->NbrObj+1][4] = 0;

      //Mise à jour de la table
      volatile int k = hsh->Head[iVer1+iVer2];
      while(hsh->LstObj[k][4]!=0){k=hsh->LstObj[k][4];}
      hsh->LstObj[k][4] = hsh->NbrObj+1;
      hsh->NbrObj++;
    }
    else // L'arête existe déjà
    {
      if(hsh->LstObj[pos][3] == 0){
        hsh->LstObj[pos][3] = iTri;
      }
    }
  } 
  else // Tête vide
  {
    hsh->Head[iVer1+iVer2] = hsh->NbrObj + 1;
    hsh->LstObj[hsh->NbrObj+1][0] = iVer1;
    hsh->LstObj[hsh->NbrObj+1][1] = iVer2;
    hsh->LstObj[hsh->NbrObj+1][2] = iTri;
    hsh->LstObj[hsh->NbrObj+1][3] = 0;
    hsh->LstObj[hsh->NbrObj+1][4] = 0;
    hsh->NbrObj++;
  }
  return 0;
}
int hash_suppr(volatile HashTable* hsh, int iVer1, int iVer2, int iTri)
{

  // to be implemented

  // ===> suppress this entry in the hash tab

  return 0;
}

int hash_count_head(volatile HashTable* hsh, int id)
{
  /*
    Fonction qui compte le nombre d'élément dans une table de hachage. 

    input  : 
        - hsh --> table de hachage.
        - id  --> pointeur vers le premier élément associé à la clé.
    output : 
        - 0 si la tête est vide, count = nombre d'élément dans la liste. 
  */
  if(id == 0){return 0;}

  int k = id, count = 0;
  while(k!= 0){count++;k = hsh->LstObj[k][4];}

  return count;
}
int* hash_biggest_head(volatile HashTable* hsh)
{
  /*
    Retrouve la tête dans laquelle il y a le plus d'éléments.

    input  : hsh --> table de hachage.
    output : array output --> output[0] = nombre maximal d'arête , output[1] = identifiant de la tête en question
  */
  int max = 0, temp = 0, id = 0;
  for(int i = 0; i < hsh->SizHead; i++)
  {
    temp = hash_count_head(hsh, hsh->Head[i]);
    if(temp == 0){continue;}
    else 
    {
      if(temp >= max)
      {
        max = temp;
        id = hsh->Head[i];
      }
    }
  } 
  
  int* output = malloc(sizeof(int)*2);
  *(output) = max; *(output+1) = id;
  
  return output;
}
double Av_colision(volatile HashTable* hsh)
{
  /*
    Fonction qui calcule le nombre moyen de colision dans la table de hachage.
    !! On exclut les tête vide !!

    input  : hsh    --> table dont on souhaite calculer le nombre moyen de colision.
    output : double --> nombre moyen de colision. 
  */
  int NEH = 0; // NEK --> Non Empty Head (nombre de tête non vide)
  int TC = 0; // TC --> Total colision (nombre total de colision)
  
  for(int i = 0; i < hsh->SizHead; i++)
  {
    if(hash_count_head(hsh, hsh->Head[i]) == 0){continue;}

    int id  = hsh->Head[i];
    TC += hash_count_head(hsh, id); // Mise à jour de TC
    NEH++; // Mise à jour de NEH
  }

  return (double)TC/(double)NEH;
}

// Méthodes sur les stacks 2d
stack_edg* stack_edg_init(int size)
{
  stack_edg* s = malloc(sizeof(stack_edg));

  s->list = malloc(sizeof(int2d)*(size+1));
  for(int i = 0; i <= size; i++){s->list[i][0] = 0; s->list[i][1] = 0;}

  s->top = 0;
  s->size = size;

  return s;
}
void stack_edg_add(stack_edg* s, int2d* new)
{
  /* Ajoute un élément à la stack. */

  s->top++;
  s->list[s->top][0] = new[0];
  s->list[s->top][1] = new[1];

  return;
}
void stack_edg_del(stack_edg* s)
{
  /* Retire la tête de la stack */

  s->list[s->top][0] = 0; 
  s->list[s->top][1] = 0;
  s->top--;

  return;
}
void stack_edg_cout(stack_edg* s)
{
  printf("\nCurrent top : %d\n", s->top);
  for(int i = 0; i <= s->top; i++){printf("\n%d : (%d,%d)\n", i, s->list[i][0], s->list[i][1]);}
  return;
}

// Méthodes sur les stacks
stack* stack_init(int size)
{
  stack* s = malloc(sizeof(stack));

  s->array = malloc(sizeof(int)*(size+1));
  for(int i = 0; i <= size; i++){s->array[i] = 0;}

  s->top = 0;
  s->size = size;

  return s;
}
void stack_add(stack* s, int new)
{
  /* Ajoute un élément à la stack. */

  s->top++;
  s->array[s->top] = new;

  return;
}
void stack_del(stack* s)
{
  /* Retire la tête de la stack */

  s->array[s->top] = 0;
  s->top--;

  return;

}
void stack_cout(stack* s)
{
  printf("\nCurrent top : %d\n", s->top);
  for(int i = 0; i <= s->top; i++){printf("%d : %d\n", i, s->array[i]);}
  return;
}

volatile HashTable* Hash_build(Mesh* Msh)
{

  /*
    Fonction qui construit la table de hachage des arêtes d'un maillage donné avec la clé de hachage iVer1+iVer2.

    Input  : 
        - Msh -> le maillage dont on souhaite la table de hachage.
    Output : 
        - hsh -> table de hachage des arêtes du maillage.
  */


  int iTri, iEdg, iVer1, iVer2;

  if (!Msh) return 0;

  if (Msh->TriVoi == NULL)
    Msh->TriVoi = calloc((Msh->NbrTri + 1), sizeof(int3d));

  // Initialiser la table de hachage
  volatile HashTable* hsh = hash_init(2*Msh->NbrVer, 3*Msh->NbrTri);

    // Construction effective de la table
  for(iTri = 1; iTri <= Msh->NbrTri; iTri++){
    for(iEdg = 0; iEdg < 3; iEdg++)
    {
      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];
      hash_add(hsh, iVer1, iVer2, iTri);
    }
  }

  return hsh;
}
int msh_neighbors(Mesh* Msh)
{
  /*
    Fonction qui calcule la liste des voisins des triangles du maillages en calculant la table de hachage des arêtes. 

    input  : 
        - Msh --> Maillage.
    output : 
        Pas important (la fonction agit comme un void sur Msh).
  */

  double start = clock();

  int iTri, iEdg, iVer1, iVer2, tri;
  volatile HashTable* hsh = Hash_build(Msh);

  // for(int k = 0; k < hsh->SizHead; k++){printf(" %d ", hsh->Head[k]);}

  //--- Compute the neighbors using the hash table
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) 
  {
    for (iEdg = 0; iEdg < 3; iEdg++) 
    {
      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]]; 
      int id = hsh->Head[iVer1+iVer2];

      while(id != 0)
      {
        if((hsh->LstObj[id][0] == iVer1 && hsh->LstObj[id][1] == iVer2) || (hsh->LstObj[id][0] == iVer2 && hsh->LstObj[id][1] == iVer1)) // Identifier l'arête 
        { 
          tri = (hsh->LstObj[id][2] == iTri) ? hsh->LstObj[id][3] : hsh->LstObj[id][2]; // Ne pas prendre le triangle sur lequel on se trouve
          
          // Vérifier que le triangle n'a pas déjà été stocké
          if(tri != iTri && tri != 0) // Si le triangle est valide, alors il est "stockable"
          {
            for(int p = 0; p < 3; p++)
            {
              if(Msh->TriVoi[iTri][p] == tri){break;} // Triangle déjà stocké
              if(Msh->TriVoi[iTri][p] == 0) // Voisin pas encore assigné
                {Msh->TriVoi[iTri][p] = tri; break;} // Stockage et sortie du for
              else{continue;} // 
            }
            id = hsh->LstObj[id][4]; // Passer au chaînon suivant

            // printf("\nTriangle %d : %d \t %d \t %d \n", iTri, Msh->TriVoi[iTri][0], Msh->TriVoi[iTri][1], Msh->TriVoi[iTri][2]);
          }
          else{id = hsh->LstObj[id][4];} // Passer directement au chaînon suivant
        }
        else{id = hsh->LstObj[id][4];}
      }
    }


  }
  double stop = clock();
  printf("Table de hachage calculée en %lg seconde(s) \n", (stop - start) / CLOCKS_PER_SEC);

  printf("Nombre d'arêtes récupérées : %d\n", hsh->NbrObj);

  return 1;
}
int Edges_build(Mesh* Msh)
{
  /*
    Fonction qui calcule les arêtes de bord du maillage. 

    input  : 
          - Msh --> maillage.
    output : 
          - Pas important (la fonction agit comme un void sur Msh).
  */
  
  volatile HashTable* hsh = Hash_build(Msh);

  int k = 0;
  int count = 0;
  for(int i = 1; i <= hsh->NbrObj; i++)
  {
    if(hsh->LstObj[i][3] == 0){Msh->EfrRef[k] = i;count++;k++;}
    Msh->NbrEdg = k;
  }

  printf("Nombre d'arête sur la frontière : %d", count);
  return 1;
}

void coloriage_magique(Mesh* Msh)
{ 
  /*
    A corriger : remplacer int* stack par stack* stack et réadapter l'algorithme en conséquence
  */
 
  int* stack = malloc(sizeof(int) * Msh->NbrTri); stack[0] = 0; // stack* stack = stack_init(Msh->NbrTri);
  int top = 0;
  int color = 0;
  int iTri, tri, voi;

  int* ar_com = malloc(sizeof(int)*2); 
  // Initialiser les couleurs (non colorié == -1)
  for(int i = 1; i <= Msh->NbrTri; i++){Msh->TriRef[i] = -1;}

  for(int i = 1; i <= Msh->NbrEfrMax; i++){printf("%d : %d \n", i, Msh->EfrRef[i]);}
  
  for(iTri = 1; iTri <= Msh->NbrTri; iTri++)
  {
    if(Msh->TriRef[iTri] == -1) // Le triangle n'est pas colorié
    {
      // Initialiser la pile : stack[0] = iTri 
      stack[1] = iTri; // tri = stack->array[stack->top];
      top = 1;

      while(top != 0)
      {
        tri = stack[top];
        // On commence par colorier le triangle en tête de la stack
        Msh->TriRef[tri] = color; 
        top--;
        //         Msh->TriRef[tri] = color; 
        //         stack_del(stack);
        
        // On parcourt les voisins de tri
        for(int k = 0; k < 3; k++)
        {
          voi = Msh->TriVoi[tri][k];
          if(voi == 0){continue;} // Triangle sur le bord du maillage

          ar_com = intersect(Msh->Tri[tri], Msh->Tri[voi]); // trouver l'arête commune entre les 2 triangles
          if(is_edge(Msh, ar_com, color) == 0)  // Si l'arête n'est pas dans un bord, ajouter le voisin à la stack
          {
            if(Msh->TriRef[voi] != -1){continue;} // Vérifier que le triangle n'est pas déjà colorié
            stack[top+1] = voi; // On l'ajoute à la stack
            top++; // On met le pointeur sur la tête de la stack à jour
            //             stack_add(stack, voi);
          }
        // Si l'arête est un bord, alors on ne l'ajoute pas à la stack et on continue (au début, top != 0)
        }
      }
      // En sortant du while, on a parcouru tout le domaine, donc on passe à la couleur suivante et on recommence
      color++;
    }
  }
  
  // On écrit directement dans l'attribut TriRef de Msh, donc pas d'objet à retourner.
  return; 
}
int* intersect(int3d Tri1, int3d Tri2)
{
  int* ar_com = malloc(sizeof(int)*2); ar_com[0] = 0; ar_com[1] = 0;
  int k = 0;
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      if(Tri1[i] == Tri2[j])
      {
        ar_com[k] = Tri1[i];
        k++;
        if(k == 2){break;}
      }
    }
  }

  return ar_com;
}
int is_edge(Mesh* M, int* id_vert, int color)
{
  int found = 0;

  for(int i = 1; i <= M->NbrEfr; i++)
  {
    if((M->Efr[i][0] == id_vert[0] && M->Efr[i][1] == id_vert[1]) || (M->Efr[i][0] == id_vert[1] && M->Efr[i][1] == id_vert[0])) // Identifier une arête de bord
    {
      if(M->EfrRef[i] != color)
      {
        found = 1;
        return found;
      }

    }
  }

  return found;
}

void Edges_vertices_cout(Mesh* Msh)
{
  printf("\n\n");
  for(int i = 0; i < Msh->NbrEdg; i++){printf("\n %d : %d\n",i,Msh->EfrRef[i]);}
  printf("\n\n");

  return;
}
void hash_cout_head(volatile HashTable* hsh, int Key)
{  
  /*
    Fonction qui affiche la liste associée à une certaine clé de la table de hachage. 

    input  :  
        - hsh --> table de hachage.
        - Key --> Premier pointeur sur la tête de la liste. 
    output : 
        - void
  */

  printf("Affichage de la liste associée à la clé %d\n", Key);
  
  int i = 1;
  int k  = Key;
  while(k != 0)
  {
    printf("%d : Ver1 : %d \tVer2 : %d \t Tri1 : %d \t Tri2 :%d \t Nxt : %d \n", i, hsh->LstObj[k][0], hsh->LstObj[k][1], hsh->LstObj[k][2], hsh->LstObj[k][3], hsh->LstObj[k][4]);
    k = hsh->LstObj[k][4];
    i++;
  }

  return;
}
void hash_cout(volatile HashTable* hsh)
{ 
  /*
    Fonction d'affichage de toute la table de hachage, tête par tête.

    input  : 
        - hsh --> table de hachage. 
    output : 
        - void
  */

  printf("Affichage de la table de hachage.\n");
  for (volatile int j = 0; j < hsh->SizHead; j++){printf(" %d ", hsh->Head[j]);}
  for(volatile int i = 0; i < hsh->SizHead; i++)
  {
    if(i == 0){printf("\n");}
    if(hsh->Head[i] == 0){continue;}

    printf("\n");
    hash_cout_head(hsh, hsh->Head[i]);
    printf("\n");
  }

  printf("Nombre d'objet dans la table : %d\n", hsh->NbrObj);
  printf("Nombre maximum d'objet dans la table : %d\n", hsh->NbrMaxObj);


  return;
}
void LstObj_cout(volatile HashTable* hsh)
{
  /*
    Fonction d'affichage de la table de hachage objet par objet. 

    input  : 
      - hsh --> table de hachage. 
    output : 
      - void
  */

  printf("\n\n");
  for(int j = 0; j <= hsh->NbrObj; j++){printf("%d : %d, %d, %d, %d, %d\n",j , hsh->LstObj[j][0]
                                                                             , hsh->LstObj[j][1]
                                                                             , hsh->LstObj[j][2]
                                                                             , hsh->LstObj[j][3]
                                                                             , hsh->LstObj[j][4]);}
  printf("\n\n");
  
  return;
}
void Head_cout(volatile HashTable* hsh)
{
  for(int k = 0; k < hsh->SizHead; k++){printf(" %d ", hsh->Head[k]);}
  return;
}
void TriVoi_cout(Mesh* Msh)
{
  for(int j = 0; j <= Msh->NbrTri; j++){printf("%d : %d, %d, %d\n",j , Msh->TriVoi[j][0]
                                                                     , Msh->TriVoi[j][1]
                                                                     , Msh->TriVoi[j][2]);}
  
  return;
}

int msh_write2dfield_Vertices(char* file, int nfield, double* field)
{
  int iVer;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  int sizfld[1];
  sizfld[0] = GmfSca;

  GmfSetKwd(fmsh, GmfSolAtVertices, nfield, 1, sizfld);

  for (iVer = 1; iVer <= nfield; iVer++)
    GmfSetLin(fmsh, GmfSolAtVertices, &field[iVer]);

  GmfCloseMesh(fmsh);

  return 1;
}
int msh_write2dfield_Triangles(char* file, int nfield, double* field)
{
  int iTri;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  int sizfld[1];
  sizfld[0] = GmfSca;

  GmfSetKwd(fmsh, GmfSolAtTriangles, nfield, 1, sizfld);

  for (iTri = 1; iTri <= nfield; iTri++)
    GmfSetLin(fmsh, GmfSolAtTriangles, &field[iTri]);

  GmfCloseMesh(fmsh);

  return 1;
}
int msh_write2dmetric(char* file, int nmetric, double3d* metric)
{
  int iVer;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  int sizfld[1];
  sizfld[0] = GmfSymMat;

  GmfSetKwd(fmsh, GmfSolAtVertices, nmetric, 1, sizfld);

  for (iVer = 1; iVer <= nmetric; iVer++)
    GmfSetLin(fmsh, GmfSolAtVertices, &metric[iVer][0], &metric[iVer][1], &metric[iVer][2]);

  GmfCloseMesh(fmsh);

  return 1;
}

// Fonctions de calcul des qualités d'un triangle (débuggée sur le 1er triangle du maillage carre_debug). Cf énoncé de la partie 1 pour voir les formules. 
double qual1(Mesh* M, int idTri)
{
  //double2d c0, c1, c2 = {0, 0};

  // Récupérer les coordonnée
  int loc_c0 = M->Tri[idTri][0];int loc_c1 = M->Tri[idTri][1];int loc_c2 = M->Tri[idTri][2];
  int3d L_loc = {loc_c0, loc_c1, loc_c2};

  double V = 0.5 * ((M->Crd[loc_c1][0]-M->Crd[loc_c0][0])*(M->Crd[loc_c2][1]-M->Crd[loc_c0][1])
                   -(M->Crd[loc_c2][0]-M->Crd[loc_c0][0])*(M->Crd[loc_c1][1]-M->Crd[loc_c0][1])); // Aire

  if(V < 0 || V == 0){printf("Triangle défectueux."); exit(-1);}

  // Calcul des longueur de chaque côté
  double l0 = pow(M->Crd[L_loc[0]][0]-M->Crd[L_loc[1]][0], 2) + pow(M->Crd[L_loc[0]][1]-M->Crd[L_loc[1]][1], 2); // 0.125
  double l1 = pow(M->Crd[L_loc[1]][0]-M->Crd[L_loc[2]][0], 2) + pow(M->Crd[L_loc[1]][1]-M->Crd[L_loc[2]][1], 2); // 0.125
  double l2 = pow(M->Crd[L_loc[2]][0]-M->Crd[L_loc[0]][0], 2) + pow(M->Crd[L_loc[2]][1]-M->Crd[L_loc[0]][1], 2); // 0.5

double alpha1 = sqrt(3)/12; // Alpha1 pour un triangle équilatéral de côté 1

return alpha1*(l0+l1+l2)/V;
}
double qual2(Mesh* M, int idTri)
{
  //double2d c0, c1, c2 = {0, 0};

  // Récupérer les coordonnée
  int loc_c0 = M->Tri[idTri][0];int loc_c1 = M->Tri[idTri][1];int loc_c2 = M->Tri[idTri][2];
  int3d L_loc = {loc_c0, loc_c1, loc_c2};

  double V = 0.5 * ((M->Crd[loc_c1][0]-M->Crd[loc_c0][0])*(M->Crd[loc_c2][1]-M->Crd[loc_c0][1])
                   -(M->Crd[loc_c2][0]-M->Crd[loc_c0][0])*(M->Crd[loc_c1][1]-M->Crd[loc_c0][1])); // Aire

  double l0 = sqrt(pow(M->Crd[L_loc[0]][0]-M->Crd[L_loc[1]][0], 2) + pow(M->Crd[L_loc[0]][1]-M->Crd[L_loc[1]][1], 2)); // 0.125
  double l1 = sqrt(pow(M->Crd[L_loc[1]][0]-M->Crd[L_loc[2]][0], 2) + pow(M->Crd[L_loc[1]][1]-M->Crd[L_loc[2]][1], 2)); // 0.125
  double l2 = sqrt(pow(M->Crd[L_loc[2]][0]-M->Crd[L_loc[0]][0], 2) + pow(M->Crd[L_loc[2]][1]-M->Crd[L_loc[0]][1], 2)); // 0.5
  
  double rho = 2*V/(l0+l1+l2); // Rayon du cercle inscrit

  // Arête la plus longue du triangle
  double h_max = (l0 <= l1) ? l1 : l0; 
  h_max = (h_max <= l2) ? l2 : h_max;
  
  double alpha2 = sqrt(3)/6; // Alpha2 pour un triangle equilatéral de côté 1

  return alpha2*h_max/rho; 
}

// Localisation d'un point dans le maillage
double* compute_BC(Mesh* Msh, int id_tri, double* P)
{
  /*Calcul des coordonnées barycentriques d'un point pour un triangle d'identifiant id_tri dans Msh un maillage.*/
 
  double* P0 = malloc(sizeof(double2d));
  double* P1 = malloc(sizeof(double2d));
  double* P2 = malloc(sizeof(double2d));
  
  // Récupérer les identifiants de chaque sommets du triangle
  int id_vert0 = Msh->Tri[id_tri][0]; 
  int id_vert1 = Msh->Tri[id_tri][1];
  int id_vert2 = Msh->Tri[id_tri][2];
  
  // Récupérer les coordonnées de chaque points
  P0[0] = Msh->Crd[id_vert0][0]; P0[1] = Msh->Crd[id_vert0][1];
  P1[0] = Msh->Crd[id_vert1][0]; P1[1] = Msh->Crd[id_vert1][1];
  P2[0] = Msh->Crd[id_vert2][0]; P2[1] = Msh->Crd[id_vert2][1];

  // Calcul effectif des coordonnées barycentrique
  double vol_tri =  0.5*((P1[0] - P0[0])*(P2[1] - P0[1]) - (P2[0] - P0[0])*(P1[1] - P0[1]));
  if(vol_tri <= 0){printf("Triangle invalide.\n"); return NULL;}

  double* Coord_bary = malloc(sizeof(double)*3);
  Coord_bary[0] = 0; Coord_bary[1] = 0; Coord_bary[2] = 0;

  Coord_bary[0] = (0.5/vol_tri)*((P1[0] - P[0])*(P2[1] - P[1]) - (P2[0] - P[0])*(P1[1] - P[1])); // |P P1 P2 |
  Coord_bary[1] = (0.5/vol_tri)*((P[0] - P0[0])*(P2[1] - P0[1]) - (P2[0] - P0[0])*(P[1] - P0[1])); // |P0 P P2|
  Coord_bary[2] = (0.5/vol_tri)*((P1[0] - P0[0])*(P[1] - P0[1]) - (P[0] - P0[0])*(P1[1] - P0[1])); // |P0 P1 P|

  free(P0);free(P1);free(P2);
  return Coord_bary;
}
int is_point_localised(Mesh* Msh, int id_tri, double* P)
{
  /* return 0 --> Le point n'est pas dans le triangle, return 1 --> Le point est dans le triangle */
  int k = 0; 

  double* Coord_bary_P = compute_BC(Msh, id_tri, P);

  while(k < 3)
  {
    if(Coord_bary_P[k] < 0){free(Coord_bary_P);return 0;}
    else{k++;}
  }

  free(Coord_bary_P);
  return 1;
}
int identify_voi(Mesh* Msh,int id_tri, int vert1, int vert2)
{
  /*Identifier le voisin du triangle id_tri du maillage Msh partageant les sommet vert1 et vert2.*/

  int voi = 0;

  for(int i = 0; i < 3; i++) // parcourir chaque voisin du triangle
  {
    if(Msh->TriVoi[id_tri][i] == 0){continue;}
    voi = Msh->TriVoi[id_tri][i];  
    if ((Msh->Tri[voi][0] == vert1 && Msh->Tri[voi][1] == vert2) || (Msh->Tri[voi][0] == vert2 && Msh->Tri[voi][1] == vert1)){return voi;}
    if ((Msh->Tri[voi][1] == vert1 && Msh->Tri[voi][2] == vert2) || (Msh->Tri[voi][1] == vert2 && Msh->Tri[voi][2] == vert1)){return voi;}
    if ((Msh->Tri[voi][2] == vert1 && Msh->Tri[voi][0] == vert2) || (Msh->Tri[voi][2] == vert2 && Msh->Tri[voi][0] == vert1)){return voi;}
  }

  printf("L'arête (%d,%d) ne semble partagée par aucun des voisins de %d, il doit y avoir une erreur...\n", vert1, vert2, id_tri);
  return -1;
}
int localising(Mesh* Msh, double* P)
{
  /*Pour l'instant, on fait l'hypothèse que le maillage est convexe.*/
  
  int itri = 1, count = 0;
  double* bary = malloc(sizeof(double)*3); for(int i = 0; i < 3; i++){bary[i] = 0;} 

  while(is_point_localised(Msh, itri, P) != 1) // Tant que le point n'est pas dans le triangle itri
  {
    // On récupère les coordonnées barycentriques
    bary = compute_BC(Msh, itri, P);

    // Les calculs conditionnels ci-dessus dépendent fortement de l'ordre dans lequel les beta0, beta1 et beta2 sont stockés
    if(bary[0] < 0)
    {
      itri = identify_voi(Msh, itri, Msh->Tri[itri][1], Msh->Tri[itri][2]); // Assigner le bon voisin
      if(itri == -1){printf("Erreur dans l'idenfication du voinsin.\n");exit(-1);}
      for(int i = 0; i < 3; i++){bary[i] = 0;} // Réinitialiser le tableau des coordonnées barycentriques par sécurité
    }
    else 
    {
      if(bary[0] >= 0 && bary[1] < 0)
      {
        itri = identify_voi(Msh,itri,Msh->Tri[itri][0], Msh->Tri[itri][2]);
        if(itri == -1){printf("Erreur dans l'idenfication du voinsin.\n");exit(-1);}
        for(int i = 0; i < 3; i++){bary[i] = 0;} 
      }
      if(bary[0] >= 0 && bary[1] >= 0)
      {
        itri = identify_voi(Msh, itri, Msh->Tri[itri][0], Msh->Tri[itri][1]);
        if(itri == -1){printf("Erreur dans l'idenfication du voinsin.\n");exit(-1);}
        for(int i = 0; i < 3; i++){bary[i] = 0;} 
      }
    }

    count++;
  }

  printf("\nPoint localisé en %d étape(s).\n", count);
  free(bary);

  return itri;
}

// Critère de Delaunay
double det_tri_points(double* P0, double* P1, double* P2)
{
  return 0.5*((P1[0] - P0[0])*(P2[1] - P0[1]) - (P2[0] - P0[0])*(P1[1] - P0[1])); // |P P1 P2 |
}
double det_tri(Mesh* Msh, int id_tri)
{
  /*Calculer l'aire signée d'un triangle d'identifiant id_tri dans le maillage Msh.*/
  
  // Récupérer les coordonnée
  int loc_c0 = Msh->Tri[id_tri][0];int loc_c1 = Msh->Tri[id_tri][1];int loc_c2 = Msh->Tri[id_tri][2];

  double V = 0.5 * ((Msh->Crd[loc_c1][0]- Msh->Crd[loc_c0][0])*(Msh->Crd[loc_c2][1] - Msh->Crd[loc_c0][1])
                   -(Msh->Crd[loc_c2][0] - Msh->Crd[loc_c0][0])*(Msh->Crd[loc_c1][1] - Msh->Crd[loc_c0][1])); // Aire

  return V;
}
double circ_circle_ray(Mesh* Msh, int id_tri)
{
  /*Calculer le rayon du cercle inscrit d'un triangle d'identifiant id_tri dans un maillage Msh.*/

  // Récupérer les coordonnée
  int loc_c0 = Msh->Tri[id_tri][0], loc_c1 = Msh->Tri[id_tri][1], loc_c2 = Msh->Tri[id_tri][2];
  int3d L_loc = {loc_c0, loc_c1, loc_c2};

  double V = det_tri(Msh, id_tri); // Aire

  // Récupérer la longeur de chaque côté du triangle
  double l0 = sqrt(pow(Msh->Crd[L_loc[0]][0] - Msh->Crd[L_loc[1]][0], 2) + pow(Msh->Crd[L_loc[0]][1] - Msh->Crd[L_loc[1]][1], 2));
  double l1 = sqrt(pow(Msh->Crd[L_loc[1]][0] - Msh->Crd[L_loc[2]][0], 2) + pow(Msh->Crd[L_loc[1]][1] - Msh->Crd[L_loc[2]][1], 2));
  double l2 = sqrt(pow(Msh->Crd[L_loc[2]][0] - Msh->Crd[L_loc[0]][0], 2) + pow(Msh->Crd[L_loc[2]][1] - Msh->Crd[L_loc[0]][1], 2));

  // Calcul effectif du rayon du cercle circonscrit
  double ray = (l0*l1*l2)/(4*V);

  return ray;
}
double* centre_circ_circle(Mesh* Msh, int id_tri)
{
  /* On calcule ici :
        | x0**2 + y0**2  x0  1 |
  x_O = | x1**2 + y1**2  X1  1 |  
        | x2**2 + y2**2  x2  1 |
  
        | x0**2 + y0**2  y0  1 |
  y_O = | x1**2 + y1**2  y1  1 |  
        | x2**2 + y2**2  y2  1 |

  et (x_O,y_O) est le centre du circle circonscrit du triangle d'identifiant id_tri dans le maillage Msh.
  On pourra retrouver la formule utilisée pour calculer x_O et y_O en développant les determinants 
  ci-dessus suivant la première colonne.
  */

  // Récupérer les coordonnées des points
  double* P0 = malloc(sizeof(double2d));
  double* P1 = malloc(sizeof(double2d));
  double* P2 = malloc(sizeof(double2d));

  int id_vert0 = Msh->Tri[id_tri][0]; 
  int id_vert1 = Msh->Tri[id_tri][1];
  int id_vert2 = Msh->Tri[id_tri][2];
  
  P0[0] = Msh->Crd[id_vert0][0]; P0[1] = Msh->Crd[id_vert0][1];
  P1[0] = Msh->Crd[id_vert1][0]; P1[1] = Msh->Crd[id_vert1][1];
  P2[0] = Msh->Crd[id_vert2][0]; P2[1] = Msh->Crd[id_vert2][1];
  
  // Calcul des valeurs de la première colonnes
  double a = pow(P0[0],2) + pow(P0[1],2), b = pow(P1[0],2) + pow(P1[1],2), c = pow(P2[0],2) + pow(P2[1],2);

  // Calcul de l'air du triangle : 0.5*det(tri)
  double v = det_tri(Msh, id_tri);
  v = (v < 0) ? -v : v; // prendre la valeur absolue

  // Calcul de x_O et y_O
  double x_O = (1/(4*v)) * (a*(P1[1] - P2[1]) - b*(P0[1] - P2[1]) + c*(P0[1] - P1[1]));
  double y_O = (-1/(4*v)) * (a*(P1[0] - P2[0]) - b*(P0[0] - P2[0]) + c*(P0[0] - P1[0]));

  double* centre = malloc(sizeof(double)*2);
  centre[0] = x_O; centre[1] = y_O;

  free(P0); free(P1); free(P2);
  return centre;
}
int empty_sphere_criterion(Mesh* Msh, int id_tri, double* P)
{
  /*
    Vérifier que le point à deux coordonnées P vérifie le critère de la sphère vide pour le triangle
    d'identifiant id_tri dans le maillage Msh.

    Pour chaque voisin du triangle id_tri, on calcul le centre et le rayon du cercle circonscrit et on vérifie que la
    distance en P et le centre est supérieur au rayon. Si c'est le cas, alors int bool == 1, sinon bool == 0.
  */

  // Init
  int bool = 0;
  double ray = circ_circle_ray(Msh, id_tri);
  double* centre = centre_circ_circle(Msh, id_tri);

  // Calcul de la distance et vérification du critère

  double dist = sqrt(pow(P[0] - centre[0],2) + pow(P[1] - centre[1], 2));

  bool = (dist < ray) ? 0 : 1;

  return bool;
}

// Insertion d'un point avec création de sa cavité suivant le critère de Delaunay et mise à jour du maillage
void push_NbrVerMax(Mesh* Msh)
{
  /*Mettre à jour le nombre de point maximum dans le maillage.*/

  Msh->NbrVerMax++;
  return;
}
int* compute_cavity(Mesh* Msh, int id_tri, double* P)
{
  // Initialiser les variables qui seront utilisées
  int* cavity = malloc(sizeof(int)*20); // Array qui stocvkera les indices des triangles qui seront dans la cavité
  for(int i = 0; i < 20; i++){cavity[i] = 0;} // Mettre le array à 0 (ce sera utile dans la fonction insertion)
  int index = 0, tri = 0, already_stored = 0; 

  // Initialiser la stack
  stack* stack = stack_init(20);
  stack->top++;
  stack->array[1] = id_tri;

  while(stack->top!=0)
  {
    tri = stack->array[stack->top]; // Stocker le triangle courant
    stack_del(stack);               // Dépiler

    // Vérifier que le triangle courant n'est pas déjà stocké
    for(int l = 0; l < index; l++)
    {
      if(cavity[l] == tri){already_stored = 1;}
    }
    
    if(empty_sphere_criterion(Msh, tri, P) == 0) // Le point est dans le cercle circonscrit du triangle id_tri
    {
      if(already_stored == 0)
      {
        cavity[index] = tri;
        index++;
        for(int k = 0; k < 3; k++)
        {
          // Empiler les voisins qui ne sont pas sur le bord
          if(Msh->TriVoi[tri][k] != 0){stack_add(stack, Msh->TriVoi[tri][k]);}
          else{continue;}
        }
      }
      else{already_stored = 0;}
    } 
  }
  for(int i = 0; i < stack->size; i++){printf("\n%d : %d\n", i, cavity[i]);}
  return cavity;
}

void insertion(Mesh* Msh, double* P)
{ 
  if(Msh->NbrVer+1 > Msh->NbrVerMax){printf("Nombre de points de maillage autorisé dépassé. Le maillage ne sera pas modifié.\n"); return;}
  
  int vert = 0, j = 0, already_stored = 0, tri = 0, ver = 0; 
  int2d* bat_cavity = NULL;

  // Localiser le point P dans le maillage
  int tri_P = localising(Msh, P);

  // Vérifier que le point n'est pas déjà dans le maillage (peut-être considérer les coordonnées barycentriques, mais je pense que ça ne fait pas beaucoup de différences en terme de perf).
  for(int i = 0; i < 3; i++)
  {
    vert = Msh->Tri[tri_P][i];
    if((Msh->Crd[vert][0] == P[0] && Msh->Crd[vert][1] == P[1]) || (Msh->Crd[vert][0] == P[1] && Msh->Crd[vert][1] == P[0]))
    {
      printf("Le point P = (%f,%f) est déjà dans le maillage. Ce dernier n'a par conséquent pas été mit à jours.\n", P[0],P[1]);
      return;
    }
  }

  // Ajouter le nouveau point au maillage
  double2d* tmp = realloc(Msh->Crd,  sizeof(double2d)*(Msh->NbrVer + 1));
  Msh->Crd = tmp;
  Msh->Crd[Msh->NbrVer+1][0] = P[0]; Msh->Crd[Msh->NbrVer+1][1] = P[1];

  // Calculer la cavité
  int* bat_cavity = compute_cavity(Msh, tri_P, P);

  Mesh* bat_mesh = msh_init();
  
  // Récupérer les noeuds de la cavité
  int* ver_cav = malloc(sizeof(int)*20);
  for(int i = 0; i < 20; i++){ver_cav[i] = 0;}
  
  while(bat_cavity[tri] != 0)
  { 
    for(int l = 0; l < 3; l++)
    {
      ver = Msh->Tri[bat_cavity[tri]][l];
      for(int m = 0; m < j; m++)
      {
        if(ver == ver_cav[m]){already_stored = 1;break;}
      }
      if(already_stored == 1){already_stored = 0;continue;} // Sommet déjà stocké
      else
      {
        ver_cav[j] = ver;
        printf("\nindex %d : %d\n", j, ver_cav[j]);
        j++;
      }
    }
    tri++;
  }
  for(int i = 0; i < j; i++)
  {
    printf("\n%d : %d\n", i, ver_cav[i]);
  }
  /* 
    tri == nombre de triangle que l'on a retiré dans le maillage, j == nombre de triangle qu'il faut ajouter
    j - tri == nombre de nouvels emplacement à créer
  */ 

  int3d* tmp_tri = realloc(Msh->Tri, sizeof(int3d)*(Msh->NbrTri + j - tri));
  Msh->Tri = tmp_tri;

  // Maintenant, on étoile pour créer les nouveaux triangles et on met à jour le maillage 


//   stack* stack = stack_init(10);
//   stack_add(stack, bat_cavity[0]);
  
//   int k = 0, cur_tri = 0;
//   while(stack->top != 0)
//   {
//     cur_tri = stack->array[stack->top];
//     stack_del(stack);

//     for(int l = 0; l < 3; l++)
//     {
//       int ver1 = Msh->Tri[cur_tri][l];

//     }
//     if(k < tri)
//     {
//       if(det_tri_points(Msh->Crd[ver_cav[k]], Msh->Crd[ver_cav[k+1]], P) < 0)  // Vérifier que le triangle est bien orienté
//       {
//         // Triangle (ar_com[0],ar_com[1],P) mal orienté, on ajoute alors le triangle (ar_com[0],P,ar_com[1])

//         Msh->Tri[bat_cavity[k]][0] = ver_cav[k];
//         Msh->Tri[bat_cavity[k]][1] = Msh->NbrVer+1;
//         Msh->Tri[bat_cavity[k]][2] = ver_cav[k+1];
//         k++;
//       }
//       else
//       {
//         // Triangle (ar_com[0],ar_com[1],P) bien orienté, on peut l'ajouter tel quel
//         Msh->Tri[bat_cavity[k]][0] = ver_cav[k];
//         Msh->Tri[bat_cavity[k]][1] = ver_cav[k+1];
//         Msh->Tri[bat_cavity[k]][2] = Msh->NbrVer+1;
//         k++;
//       }
//     }
//     else 
//     {
//       if(det_tri_points(Msh->Crd[ver_cav[k]], Msh->Crd[ver_cav[k]], P) < 0) 
//       {
//         // Triangle (ar_com[0],ar_com[1],P) mal orienté, on ajoute alors le triangle (ar_com[0],P,ar_com[1])
//         Msh->Tri[Msh->NbrTri+1][0] = ver_cav[k];
//         Msh->Tri[Msh->NbrTri+1][1] = Msh->NbrVer+1;
//         Msh->Tri[Msh->NbrTri+1][2] = ver_cav[k+1];
        
//         // Mettre à jour le nombre de triangle dans le maillage
//         Msh->NbrTri++;
//         k++;
//       }
//       else
//       {
//         // Triangle (ar_com[0],ar_com[1],P) bien orienté, on peut l'ajouter tel quel
//         Msh->Tri[Msh->NbrTri+1][0] = ver_cav[k];
//         Msh->Tri[Msh->NbrTri+1][1] = ver_cav[k+1];
//         Msh->Tri[Msh->NbrTri+1][2] = Msh->NbrVer+1;

//         // Mettre à jour le nombre de triangle dans le maillage
//         Msh->NbrTri++;
//         k++;
//       }
//     }
//   }

//   // Mise à  jour final : on incrémente le nombre de noeud dans le maillage
//   Msh->NbrVer++;
//   printf("\nLe point P = (%f,%f) à été inséré avec succès !\n", P[0], P[1]);
    
//   return;
// }
