#include "cylc_menu.h"
#include "cylc_local_menu.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <complex.h>
#define PI 3.14159265358979323846

double *transit_flux_drop_claret_cython(double *flux, double a1, double a2, double a3, double a4, double rp_over_rs, double *z, int size){

    if (size == 0){return z;}

    double *G1 = malloc(30 * sizeof(double));
    double *G2 = malloc(30 * sizeof(double));

    G1[0] = 0.1028526528935588;
    G1[1] = 0.1028526528935588;
    G1[2] = 0.1017623897484055;
    G1[3] = 0.1017623897484055;
    G1[4] = 0.0995934205867953;
    G1[5] = 0.0995934205867953;
    G1[6] = 0.0963687371746443;
    G1[7] = 0.0963687371746443;
    G1[8] = 0.0921225222377861;
    G1[9] = 0.0921225222377861;
    G1[10] = 0.0868997872010830;
    G1[11] = 0.0868997872010830;
    G1[12] = 0.0807558952294202;
    G1[13] = 0.0807558952294202;
    G1[14] = 0.0737559747377052;
    G1[15] = 0.0737559747377052;
    G1[16] = 0.0659742298821805;
    G1[17] = 0.0659742298821805;
    G1[18] = 0.0574931562176191;
    G1[19] = 0.0574931562176191;
    G1[20] = 0.0484026728305941;
    G1[21] = 0.0484026728305941;
    G1[22] = 0.0387991925696271;
    G1[23] = 0.0387991925696271;
    G1[24] = 0.0287847078833234;
    G1[25] = 0.0287847078833234;
    G1[26] = 0.0184664683110910;
    G1[27] = 0.0184664683110910;
    G1[28] = 0.0079681924961666;
    G1[29] = 0.0079681924961666;

    G2[0] = -0.0514718425553177;
    G2[1] = 0.0514718425553177;
    G2[2] = -0.1538699136085835;
    G2[3] = 0.1538699136085835;
    G2[4] = -0.2546369261678899;
    G2[5] = 0.2546369261678899;
    G2[6] = -0.3527047255308781;
    G2[7] = 0.3527047255308781;
    G2[8] = -0.4470337695380892;
    G2[9] = 0.4470337695380892;
    G2[10] = -0.5366241481420199;
    G2[11] = 0.5366241481420199;
    G2[12] = -0.6205261829892429;
    G2[13] = 0.6205261829892429;
    G2[14] = -0.6978504947933158;
    G2[15] = 0.6978504947933158;
    G2[16] = -0.7677774321048262;
    G2[17] = 0.7677774321048262;
    G2[18] = -0.8295657623827684;
    G2[19] = 0.8295657623827684;
    G2[20] = -0.8825605357920527;
    G2[21] = 0.8825605357920527;
    G2[22] = -0.9262000474292743;
    G2[23] = 0.9262000474292743;
    G2[24] = -0.9600218649683075;
    G2[25] = 0.9600218649683075;
    G2[26] = -0.9836681232797472;
    G2[27] = 0.9836681232797472;
    G2[28] = -0.9968934840746495;
    G2[29] = 0.9968934840746495;

    double A1 = - (2.0 * (1.0 - a1 - a2 - a3 - a4) / 4);
    double A2 = - (2.0 * a1 / 5);
    double A3 = - (2.0 * a2 / 6);
    double A4 = - (2.0 * a3 / 7);
    double A5 = - (2.0 * a4 / 8);
    double A15 = A1 + A2 + A3 + A4 + A5;

    int ii;
    double z_over_rs;
    double zsq;
    double sum_z_rp;
    double rpsq;
    double dif_z_rprs;
    double sqr_dif_z_rprs;
    double plus;
    double minus;
    double star;
    double test;
    double ph;
    double th;
    double total_flux;
    double r1mu44;
    double r1mu24;
    double r1mu14;
    double integral_r_rprs;
    double integral_r_0;

    rpsq = rp_over_rs * rp_over_rs;
    r1mu44 = 1.0 - rpsq;
    r1mu24 = sqrt(r1mu44);
    r1mu14 = sqrt(r1mu24);
    integral_r_rprs = (A1 * r1mu44 + A2 * r1mu44 * r1mu14 + A3 * r1mu44 * r1mu24 + A4 * r1mu44 * r1mu24 * r1mu14 + A5 * r1mu44 * r1mu44);
    integral_r_0 = A15;

    total_flux = - integral_r_0 * 2.0 * PI;

    for (ii=0;ii<size;ii++){

        plus = 0;
        minus = 0;
        star = 0;
        z_over_rs = z[ii];

        if (z_over_rs > 0){
            dif_z_rprs = rp_over_rs - z_over_rs;

            if (dif_z_rprs <= -1){flux[ii] = 1;}

            else {
                sum_z_rp = z_over_rs + rp_over_rs;

                if(z_over_rs < rp_over_rs && sum_z_rp > 1 && dif_z_rprs >= 1){flux[ii] = 0;}

                else {
                    zsq = z_over_rs * z_over_rs;
                    sqr_dif_z_rprs = zsq - rpsq;

                    test = (1.0 - rpsq + zsq) / (2.0 * z_over_rs);
                    if (test < -1){test=-1;} else if (test > 1){test=1;}
//                    ph_cos = test;
//                    ph_sin = sqrt(1 - ph_cos * ph_cos);
                    ph = acos(test);

                    test = rp_over_rs / z_over_rs;
                    if (test > 1) {test=1;}
//                    th_sin = test;
//                    th_cos = sqrt(1 - th_sin * th_sin);
                    th = asin(test);

                    if (z_over_rs < rp_over_rs){

                        if(sum_z_rp <= 1){
                        plus = integral_plus_core_claret(a1, a2, a3, a4, rp_over_rs, z_over_rs, 0.0, PI, G1, G2);
                        }

                        else if(dif_z_rprs < 1){
                        plus = integral_plus_core_claret(a1, a2, a3, a4, rp_over_rs, z_over_rs, ph, PI, G1, G2);
                        star = - integral_r_0 * ph;
                        }

                    }

                    else if (z_over_rs == rp_over_rs){

                        if(sum_z_rp <= 1){
                        plus = integral_plus_core_claret(a1, a2, a3, a4, rp_over_rs, z_over_rs, 0.0, 0.5 * PI, G1, G2);
                        }

                        else{
                        plus = integral_plus_core_claret(a1, a2, a3, a4, rp_over_rs, z_over_rs, ph, 0.5 * PI, G1, G2);
                        star = - integral_r_0 * ph;
                        }

                    }

                    else{

                        if(sum_z_rp <= 1){
                        plus = integral_plus_core_claret(a1, a2, a3, a4, rp_over_rs, z_over_rs, 0.0, th, G1, G2);
                        minus = integral_minus_core_claret(a1, a2, a3, a4, rp_over_rs, z_over_rs, 0.0, th, G1, G2);
                        }

                        else{

                            if(sqr_dif_z_rprs < 1){
                            plus = integral_plus_core_claret(a1, a2, a3, a4, rp_over_rs, z_over_rs, ph, th, G1, G2);
                            minus = integral_minus_core_claret(a1, a2, a3, a4, rp_over_rs, z_over_rs, 0.0, th, G1, G2);
                            star = - integral_r_0 * ph;
                            }

                            else if(sqr_dif_z_rprs == 1){
                            minus = integral_minus_core_claret(a1, a2, a3, a4, rp_over_rs, z_over_rs, 0.0, th, G1, G2);
                            star = - integral_r_0 * ph;
                            }

                            else if(-1 < dif_z_rprs){
                            minus = integral_minus_core_claret(a1, a2, a3, a4, rp_over_rs, z_over_rs, 0.0, ph, G1, G2);
                            star = - integral_r_0 * ph;
                            }

                        }

                    }

                    flux[ii] = 1 - (2.0 / total_flux) * (plus + star - minus);

                }

            }

        }

        else if (z_over_rs < 0){flux[ii] = 1;}

        else {
            if (rp_over_rs < 1){flux[ii] = 1 + (integral_r_rprs - integral_r_0) / integral_r_0;}
            else {flux[ii] = 0;}
        }

    }

return ;
}

