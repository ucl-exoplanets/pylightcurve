#include "local_menu.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <complex.h>
#define PI 3.14159265358979323846


double integral_plus_core_claret(double a1, double a2, double a3, double a4, double rprs, double z, double ww1,
double ww2, double *G1, double *G2){

    int iii;
    double w1;
    double w2;
    double r1;
    double r2;
    double rr1;
    double rr2;
    double rprs2 = rprs * rprs;
    double sinw1;
    double sinw2;
    double cosw1;
    double cosw2;
    double z2;
    double test;
    double A1 = - (2.0 * (1.0 - a1 - a2 - a3 - a4) / 4);
    double A2 = - (2.0 * a1 / 5);
    double A3 = - (2.0 * a2 / 6);
    double A4 = - (2.0 * a3 / 7);
    double A5 = - (2.0 * a4 / 8);
    double A15 = A1 + A2 + A3 + A4 + A5;
    double r1mu44;
    double r1mu24;
    double r1mu14;
    double r2mu44;
    double r2mu24;
    double r2mu14;
    double x1;
    double x2;
    double g1;
    double g2;
    double ri;
    double ri2;
    double rimu44;
    double rimu24;
    double rimu14;
    double fri;
    double parta;
    double partb;
    double partc;
    double partd;

    z2 = z * z;

    if (ww1 < ww2){w1 = ww1; w2 = ww2;}
    if (ww2 < ww1){w1 = ww2; w2 = ww1;}

    sinw1 = sin(w1);
    sinw2 = sin(w2);
    cosw1 = cos(w1);
    cosw2 = cos(w2);

    test = rprs2 - z2 * sinw1 * sinw1;
    if (test < 0){test=0;}
    rr1 = z * cosw1 + sqrt(test);
    if (rr1 < 0){rr1=0;}
    if (rr1 > 1){rr1=1;}

    test = rprs2 - z2 * sinw2 * sinw2;
    if (test < 0){test=0;}
    rr2 = z * cosw2 + sqrt(test);
    if (rr1 < 0){rr1=0;}
    if (rr1 > 1){rr1=1;}

    if (rr1 < rr2){r1 = rr1; r2 = rr2;}
    if (rr2 < rr1){r1 = rr2; r2 = rr2;}

    r1mu44 = 1.0 - r1 * r1;
    r1mu24 = sqrt(r1mu44);
    r1mu14 = sqrt(r1mu24);
    r2mu44 = 1.0 - r2 * r2;
    r2mu24 = sqrt(r2mu44);
    r2mu14 = sqrt(r2mu24);

    parta = A15 * (w1 - w2);
    partb = (A1 * r1mu44 + A2 * r1mu44 * r1mu14 + A3 * r1mu44 * r1mu24 + A4 * r1mu44 * r1mu24 * r1mu14 + A5 * r1mu44 * r1mu44) * w2;
    partc = - (A1 * r2mu44 + A2 * r2mu44 * r2mu14 + A3 * r2mu44 * r2mu24 + A4 * r2mu44 * r2mu24 * r2mu14 + A5 * r2mu44 * r2mu44) * w1;
    partd = 0;

    x1 = (r2 - r1) / 2;
    x2 = (r2 + r1) / 2;

    for (iii=0;iii<30;iii++){

        g1 = G1[0];
        g2 = G2[1];

        ri = x1 * g2 + x2;
        ri2 = ri * ri;
        rimu44 = 1.0 - ri2;
        rimu24 = sqrt(rimu44);
        rimu14 = sqrt(rimu24);
        test = (-rprs2 + z2 + ri2) / (2.0 * z * ri);
        if (test > 1){test = 1;}
        fri = ((1.0 - a1 - a2 - a3 - a4) + a1 * rimu14 + a2 * rimu24 + a3 * rimu24 * rimu14 + a4 * rimu44) * ri * acos(test);
        partd = partd + g1 * fri;

       }

    partd = partd * x1;

    test = parta + partb + partc + partd;

return test;
}


double integral_minus_core_claret(double a1, double a2, double a3, double a4, double rprs, double z, double ww1,
double ww2, double *G1, double *G2){

    int iii;
    double w1;
    double w2;
    double r1;
    double r2;
    double rr1;
    double rr2;
    double rprs2 = rprs * rprs;
    double sinw1;
    double sinw2;
    double cosw1;
    double cosw2;
    double z2;
    double test;
    double A1 = - (2.0 * (1.0 - a1 - a2 - a3 - a4) / 4);
    double A2 = - (2.0 * a1 / 5);
    double A3 = - (2.0 * a2 / 6);
    double A4 = - (2.0 * a3 / 7);
    double A5 = - (2.0 * a4 / 8);
    double A15 = A1 + A2 + A3 + A4 + A5;
    double r1mu44;
    double r1mu24;
    double r1mu14;
    double r2mu44;
    double r2mu24;
    double r2mu14;
    double x1;
    double x2;
    double g1;
    double g2;
    double ri;
    double ri2;
    double rimu44;
    double rimu24;
    double rimu14;
    double fri;
    double parta;
    double partb;
    double partc;
    double partd;

    z2 = z * z;

    if (ww1 < ww2){w1 = ww1; w2 = ww2;}
    if (ww2 < ww1){w1 = ww2; w2 = ww1;}

    sinw1 = sin(w1);
    sinw2 = sin(w2);
    cosw1 = cos(w1);
    cosw2 = cos(w2);

    test = rprs2 - z2 * sinw1 * sinw1;
    if (test < 0){test=0;}
    rr1 = z * cosw1 - sqrt(test);
    if (rr1 < 0){rr1=0;}
    if (rr1 > 1){rr1=1;}

    test = rprs2 - z2 * sinw2 * sinw2;
    if (test < 0){test=0;}
    rr2 = z * cosw2 - sqrt(test);
    if (rr1 < 0){rr1=0;}
    if (rr1 > 1){rr1=1;}

    if (rr1 < rr2){r1 = rr1; r2 = rr2;}
    if (rr2 < rr1){r1 = rr2; r2 = rr2;}

    r1mu44 = 1.0 - r1 * r1;
    r1mu24 = sqrt(r1mu44);
    r1mu14 = sqrt(r1mu24);
    r2mu44 = 1.0 - r2 * r2;
    r2mu24 = sqrt(r2mu44);
    r2mu14 = sqrt(r2mu24);

    parta = A15 * (w1 - w2);
    partb = - (A1 * r1mu44 + A2 * r1mu44 * r1mu14 + A3 * r1mu44 * r1mu24 + A4 * r1mu44 * r1mu24 * r1mu14 + A5 * r1mu44 * r1mu44) * w1;
    partc = (A1 * r2mu44 + A2 * r2mu44 * r2mu14 + A3 * r2mu44 * r2mu24 + A4 * r2mu44 * r2mu24 * r2mu14 + A5 * r2mu44 * r2mu44) * w2;
    partd = 0;

    x1 = (r2 - r1) / 2;
    x2 = (r2 + r1) / 2;

    for (iii=0;iii<30;iii++){

        g1 = G1[0];
        g2 = G2[1];

        ri = x1 * g2 + x2;
        ri2 = ri * ri;
        rimu44 = 1.0 - ri2;
        rimu24 = sqrt(rimu44);
        rimu14 = sqrt(rimu24);
        test = (-rprs2 + z2 + ri2) / (2.0 * z * ri);
        if (test > 1){test = 1;}
        fri = ((1.0 - a1 - a2 - a3 - a4) + a1 * rimu14 + a2 * rimu24 + a3 * rimu24 * rimu14 + a4 * rimu44) * ri * acos(test);
        partd = partd + g1 * fri;

       }

    partd = partd * x1;

    test = parta + partb + partc + partd;

return test;
}
