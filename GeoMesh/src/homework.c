#include "fem.h"


#include <math.h>

double distance(double x, double y, double xc, double yc) {
    return sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
}

double hermiteInterpolation(double d, double h0, double h_star, double d_star) {
    double t = d / d_star;
    if (t > 1) return h_star; // Si on est au-delà de la zone de transition
    return h0 * (1 - 3*t*t + 2*t*t*t) + h_star * (3*t*t - 2*t*t*t);
}


double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();


    double h = theGeometry->h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
  
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;


//
//     A modifier !
//     
// Your contribution starts here ....
//

// Distances aux centres des cercles
    double d_encoche = fmax(0, distance(x, y, x0, y0) - r0);
    double d_trou = fmax(0, distance(x, y, x1, y1) - r1);


    // Interpolation d'Hermite
    double h_encoche = hermiteInterpolation(d_encoche, h0, h, d0);
    double h_trou = hermiteInterpolation(d_trou, h1, h, d1);

    // Retourne la plus petite valeur (raffinement optimal)
    return fmin(h_encoche, h_trou);
    
    
//   
// Your contribution ends here :-)
//

}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate; //Lx
    double h = theGeometry->LyPlate; //Ly
     
    double x0 = theGeometry->xNotch; // xNotch
    double y0 = theGeometry->yNotch; // yNotch
    double r0 = theGeometry->rNotch; // rNotch
    
    
    double x1 = theGeometry->xHole; //xHole
    double y1 = theGeometry->yHole; //yHole
    double r1 = theGeometry->rHole; //rHole
 
//
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
//      On soustrait les cercles du rectangle :-)
//
 
    int ierr;

    int idPlate = gmshModelOccAddRectangle(-(w/2), -(h/2), 0.0, w, h, 0, -1, &ierr);   
    printf("idPlate = %d\n", idPlate);
    ErrorGmsh(ierr);

    int idNotch = gmshModelOccAddDisk(x0, y0, 0.0, r0, r0, 1, NULL, 0, NULL, 0, &ierr); 
    printf("idNotch = %d\n", idNotch);
    ErrorGmsh(ierr);

    int idHole = gmshModelOccAddDisk(x1, y1, 0.0, r1, r1, 2, NULL, 0, NULL, 0, &ierr);    
    printf("idHole = %d\n", idHole);
    ErrorGmsh(ierr);


    
    int plate[] = {2,idPlate};// 2 signifie surface ==> défini la plaque commme une surface
    int notch[] = {2,idNotch};
    int hole[]  = {2,idHole};
    gmshModelOccCut(plate,2,notch,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    ErrorGmsh(ierr);
    gmshModelOccCut(plate,2,hole,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    ErrorGmsh(ierr);


 
//
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr);
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);
//    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
//   gmshFltkInitialize(&ierr);
//   gmshFltkRun(&ierr);  chk(ierr);
//
    
}