#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmsh/gmsh-4.13.1-Linux64-sdk/include/gmsh.h>

#define N 50  // Nombre de points pour le profil

void generate_naca2412(double x[], double y[], int n) {
    double m = 0.02, p = 0.4, t = 0.12;
    
    for (int i = 0; i < n; i++) {
        x[i] = (double)i / (n - 1);
        double yt = 5 * t * (0.2969 * sqrt(x[i]) - 0.1260 * x[i] - 0.3516 * x[i] * x[i] +
                             0.2843 * pow(x[i], 3) - 0.1015 * pow(x[i], 4));
        double yc, theta;

        if (x[i] < p)
            yc = (m / (p * p)) * (2 * p * x[i] - x[i] * x[i]);
        else
            yc = (m / ((1 - p) * (1 - p))) * ((1 - 2 * p) + 2 * p * x[i] - x[i] * x[i]);

        theta = atan((x[i] < p) ? (2 * m / (p * p)) * (p - x[i]) : (2 * m / ((1 - p) * (1 - p))) * (p - x[i]));

        y[i] = yc + yt * cos(theta);  // Extrados
        y[n + i] = yc - yt * cos(theta);  // Intrados
    }
}

int main() {
    gmshInitialize(0, NULL, 0, NULL);

    gmshModelAdd("NACA2412");

    double x[2 * N], y[2 * N];
    generate_naca2412(x, y, N);

    int pointTags[2 * N];

    for (int i = 0; i < 2 * N; i++) {
        pointTags[i] = gmshModelGeoAddPoint(x[i], y[i], 0, 1.0, i + 1);
    }

    int curveTags[2 * N - 1];
    for (int i = 0; i < 2 * N - 1; i++) {
        curveTags[i] = gmshModelGeoAddLine(pointTags[i], pointTags[i + 1], i + 1);
    }

    gmshModelGeoSynchronize();
    gmshWrite("naca2412.geo");

    gmshFltkRun();
    gmshFinalize();
    return 0;
}
