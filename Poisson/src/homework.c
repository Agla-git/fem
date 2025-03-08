#include "fem.h"

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femGeo* theGeometry = geoMeshCreate(filename);
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->geo  = theGeometry;
    femMesh *theMesh = theGeometry->theElements;
    if (theMesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theMesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theMesh->nodes->nNodes);
    return theProblem;
}

# endif
# ifndef NOPOISSONBOUNDARY






void femPoissonFindBoundaryNodes(femPoissonProblem *theProblem)
{
    femGeo* theGeometry = theProblem->geo;  
    femMesh* theEdges = theGeometry->theEdges;  
    int nBoundary = 0;  
    
    // Création d'un tableau temporaire pour stocker les nœuds de frontière
    int *boundaryNodes = malloc(theEdges->nElem * theEdges->nLocalNode * sizeof(int));
    int *isBoundaryNode = calloc(theGeometry->theElements->nodes->nNodes, sizeof(int));

    // Parcourir toutes les arêtes du maillage
    for (int i = 0; i < theEdges->nElem; i++) {
        for (int j = 0; j < theEdges->nLocalNode; j++) {
            int node = theEdges->elem[i * theEdges->nLocalNode + j]; 

            // Vérifier si ce nœud a déjà été ajouté à la liste des frontières
            if (!isBoundaryNode[node]) {
                boundaryNodes[nBoundary++] = node;
                isBoundaryNode[node] = 1;  // Marquer comme ajouté
            }
        }
    }

    // Création du domaine de la frontière (Boundary)
    femDomain *theBoundary = malloc(sizeof(femDomain));
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains, theGeometry->nDomains * sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains - 1] = theBoundary;
    
    theBoundary->nElem = nBoundary;
    theBoundary->elem = malloc(nBoundary * sizeof(int));
    for (int i = 0; i < nBoundary; i++) {
        theBoundary->elem[i] = boundaryNodes[i]; // Stocker les indices des nœuds de frontière
    }
    
    theBoundary->mesh = NULL;
    sprintf(theBoundary->name, "Boundary");

    // Libération des tableaux temporaires
    free(boundaryNodes);
    free(isBoundaryNode);
}




/*
void femPoissonFindBoundaryNodes(femPoissonProblem *theProblem)
{
    femGeo* theGeometry = theProblem->geo;  
    femMesh* theEdges = theGeometry->theEdges; 
    int nBoundary = 0;
    
    //  A completer :-)


    femDomain *theBoundary = malloc(sizeof(femDomain));
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains,theGeometry->nDomains*sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains-1] = theBoundary;
    theBoundary->nElem = nBoundary;
    theBoundary->elem = malloc(nBoundary*sizeof(int));
    theBoundary->mesh = NULL;
    sprintf(theBoundary->name,"Boundary");
 
    // A completer :-)


*/
    
# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    if (theProblem == NULL) return;

    // Libération du système linéaire
    femFullSystemFree(theProblem->system);

    // Libération de la discrétisation et de la règle d’intégration
    femDiscreteFree(theProblem->space);
    femIntegrationFree(theProblem->rule);

    // Libération de la géométrie et du maillage
    geoMeshFree(theProblem->geo);

    // Libération de la structure du problème
    free(theProblem);
}

/*

void femPoissonFree(femPoissonProblem *theProblem)
{
    if (theProblem == NULL) return;

    // Libération du système linéaire
    femFullSystemFree(theProblem->system);

    // Libération de la discrétisation et de la règle d’intégration
    femDiscreteFree(theProblem->space);
    femIntegrationFree(theProblem->rule);

    // Libération des domaines de frontière s'ils existent
    if (theProblem->geo->nDomains > 0) {
        for (int i = 0; i < theProblem->geo->nDomains; i++) {
            free(theProblem->geo->theDomains[i]->elem); // Libère le tableau des éléments
            free(theProblem->geo->theDomains[i]);       // Libère la structure du domaine
        }
        free(theProblem->geo->theDomains); // Libère le tableau des domaines
    }

    // Libération de la géométrie et du maillage
    geoMeshFree(theProblem->geo);

    // Libération de la structure du problème
    free(theProblem);
}
*/
# endif
# ifndef NOPOISSONLOCAL

/*
void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->geo->theElements;
    
    //  A completer :-)

}
*/


void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->geo->theElements;
    
    // Copier les indices des sommets de l'élément iElem
    for (int i = 0; i < theMesh->nLocalNode; i++) {
        map[i] = theMesh->elem[iElem * theMesh->nLocalNode + i]; 
    }
    
    // Copier les coordonnées des sommets
    for (int i = 0; i < theMesh->nLocalNode; i++) {
        x[i] = theMesh->nodes->X[map[i]];
        y[i] = theMesh->nodes->Y[map[i]];
    }
}


# endif
# ifndef NOPOISSONSOLVE

/*

void femPoissonSolve(femPoissonProblem *theProblem)
{
    femMesh *theMesh = theProblem->geo->theElements;
    femDomain *theBoundary = geoGetDomain(theProblem->geo, "Boundary");
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;

    if (theSpace->n > 4) Error("Unexpected discrete space size!");

    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int map[4];  // Indices globaux des nœuds d’un élément
    int nLocal = theMesh->nLocalNode;
    double jacobian, dxdxi, dxdeta, dydxi, dydeta;
    
    // Boucle sur les éléments pour l'assemblage
    for (int iElem = 0; iElem < theMesh->nElem; iElem++) {
        femPoissonLocal(theProblem, iElem, map, x, y);  // Récupère coordonnées et indices des nœuds

        double Ke[nLocal][nLocal];  // Matrice de rigidité locale
        double Fe[nLocal];  // Second membre local
        
        // Initialisation des matrices locales
        for (int i = 0; i < nLocal; i++) {
            Fe[i] = 0.0;
            for (int j = 0; j < nLocal; j++) {
                Ke[i][j] = 0.0;
            }
        }

        // Boucle sur les points d'intégration
        for (int iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            // Calcul du Jacobien et de son inverse
            dxdxi = dydeta = 0.0;
            dxdeta = dydxi = 0.0;
            for (int i = 0; i < nLocal; i++) {
                dxdxi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }
            jacobian = dxdxi * dydeta - dxdeta * dydxi;

            // Calcul des dérivées des fonctions de forme en x et y
            for (int i = 0; i < nLocal; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxi) / jacobian;
                dphidy[i] = (dphideta[i] * dxdxi - dphidxsi[i] * dxdeta) / jacobian;
            }

            // Assemblage de la matrice locale Ke et du second membre Fe
            for (int i = 0; i < nLocal; i++) {
                Fe[i] += weight * phi[i] * jacobian;  // Second membre avec f=1 (source uniforme)
                for (int j = 0; j < nLocal; j++) {
                    Ke[i][j] += weight * (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jacobian;
                }
            }
        }

        // Assemblage dans le système global
        for (int i = 0; i < nLocal; i++) {
            int I = map[i];
            femFullSystemAdd(theSystem, I, I, Ke[i][i]);
            femFullSystemAddRightHandSide(theSystem, I, Fe[i]);
            for (int j = i + 1; j < nLocal; j++) {
                int J = map[j];
                femFullSystemAdd(theSystem, I, J, Ke[i][j]);
                femFullSystemAdd(theSystem, J, I, Ke[j][i]);  // Symétrie
            }
        }
    }

    // Imposition des conditions de Dirichlet (ex: u=0 sur la frontière)
    for (int i = 0; i < theBoundary->nElem; i++) {
        int node = theBoundary->elem[i];
        femFullSystemConstrain(theSystem, node, 0.0);
    }

    // Résolution du système linéaire
    femFullSystemSolve(theSystem);
}
*/

void femPoissonSolve(femPoissonProblem *theProblem)
{
    femMesh *theMesh = theProblem->geo->theElements;
    femDomain *theBoundary = geoGetDomain(theProblem->geo, "Boundary");
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;

    if (theSpace->n > 4) Error("Unexpected discrete space size !");  

    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int iElem, iInteg, i, j, map[4];
    int nLocal = theMesh->nLocalNode;
    int nGauss = theRule->n;

    // Assemblage du système global
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femPoissonLocal(theProblem, iElem, map, x, y);

        double A_local[4][4] = {{0}}, B_local[4] = {0};

        for (iInteg = 0; iInteg < nGauss; iInteg++) {
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            double jacobian[2][2] = {{0}}, detJ;
            for (i = 0; i < nLocal; i++) {
                jacobian[0][0] += dphidxsi[i] * x[i];
                jacobian[0][1] += dphidxsi[i] * y[i];
                jacobian[1][0] += dphideta[i] * x[i];
                jacobian[1][1] += dphideta[i] * y[i];
            }
            detJ = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];

            double invJ[2][2] = {
                {jacobian[1][1] / detJ, -jacobian[0][1] / detJ},
                {-jacobian[1][0] / detJ, jacobian[0][0] / detJ}
            };

            for (i = 0; i < nLocal; i++) {
                dphidx[i] = invJ[0][0] * dphidxsi[i] + invJ[0][1] * dphideta[i];
                dphidy[i] = invJ[1][0] * dphidxsi[i] + invJ[1][1] * dphideta[i];
            }

            double dArea = fabs(detJ) * weight;

            for (i = 0; i < nLocal; i++) {
                for (j = 0; j < nLocal; j++) {
                    A_local[i][j] += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * dArea;
                }
                B_local[i] += phi[i] * dArea;  // Terme source f(x,y) pris constant à 1
            }
        }

        for (i = 0; i < nLocal; i++) {
            for (j = 0; j < nLocal; j++) {
                theSystem->A[map[i]][map[j]] += A_local[i][j];
            }
            theSystem->B[map[i]] += B_local[i];
        }
    }

    // Imposition des conditions de Dirichlet
    for (i = 0; i < theBoundary->nElem; i++) {
        int node = theBoundary->elem[i];
        femFullSystemConstrain(theSystem, node, 0.0);
    }

    // Résolution du système
    femFullSystemEliminate(theSystem);
}


# endif



