#include <stdio.h>
#include <math.h>
#include "glfem.h"

// double u_interpolate(double x[3], double eta, double xi){
//     return x[0]*(1-xi-eta) + x[1]*eta + x[2]*(1-xi-eta);
// }

// double interpolate_array(double u[3] , double xi [3] , double eta [3] , double *v){
//     for(int i = 0; i<3; i++){
//         v[i] = u_interpolate(u,eta[i],xi[i]);
//     }
//     return 0;
// }



double jaco(double x[3],double y[3]){
    return fabs((x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]));
}

double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];
    const double weight[3] = {1.0/6.0, 1.0/6.0, 1.0/6.0};
    const double xi[3] = {1.0/6.0, 1.0/6.0, 2.0/3.0};
    const double eta[3] = {1.0/6.0, 2.0/3.0, 1.0/6.0};     
     
    double jacob = jaco(x, y);

    //////////////////////////
    
    for(int i = 0;i<3;i++){
        xLoc[i] = x[0]*(1-xi[i]-eta[i]) + x[1]*xi[i] + x[2]*eta[i];
        yLoc[i] = y[0]*(1-xi[i]-eta[i]) + y[1]*xi[i] + y[2]*eta[i];
        I+= f(xLoc[i],yLoc[i])*weight[i];
    }

    
    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    return I*jacob;
}




double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
    if(n == 0){
        return integrate(x,y,f);
    }

    double new_x0= (x[0]+x[1])/2;
    double new_x1= (x[1]+x[2])/2;
    double new_x2= (x[2]+x[0])/2;

    double new_y0= (y[0]+y[1])/2;
    double new_y1= (y[1]+y[2])/2;
    double new_y2= (y[2]+y[0])/2;

    double x1[3] = {x[0],new_x0,new_x2};
    double y1[3] = {y[0],new_y0,new_y2};

    double x2[3] = {new_x0,x[1],new_x1};
    double y2[3] = {new_y0,y[1],new_y1};

    double x3[3] = {new_x2,new_x1,x[2]};
    double y3[3] = {new_y2,new_y1,y[2]};

    double x4[3] = {new_x0,new_x1,new_x2};
    double y4[3] = {new_y0,new_y1,new_y2};

    double I = integrateRecursive(x1,y1,f,n-1) + integrateRecursive(x2,y2,f,n-1) + integrateRecursive(x3,y3,f,n-1) + integrateRecursive(x4,y4,f,n-1);
    return I;
}
