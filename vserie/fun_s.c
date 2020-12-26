/**********************************************************************************************************
 *                              AC - OpenMP -- SERIE                                                      *
 *                  Rutinas que se utilizan en el modulo gengrupos_s.c                                    *
 *                                   fun_s.c                                                              *
 **********************************************************************************************************/
#include <math.h>
#include <float.h>

#include "defineg.h"    // Definiciones


/**********************************************************************************************************
 * 1 - Funcion para calcular la distancia euclidea entre dos elementos (distancia euclidea)                 *
 *          Entrada:  2 elementos con NCAR caracteristicas (por referencia)                               *
 *          Salida:   distancia (double)                                                                  *
 **********************************************************************************************************/
double gendist (float *elem1, float *elem2) {
    double acum = 0;
    for (int i = 0; i < NCAR; i++) {
        double res = elem1[i] - elem2[i];
        acum += pow(res, 2);
    }
    return sqrt(acum);
}
/**********************************************************************************************************
 * 2 - Funcion para calcular el grupo (cluster) mas cercano (centroide mas cercano)                       *
 *          Entrada:  nelem  numero de elementos, int                                                     *
 *                    elem   elementos, una matriz de tamano MAXE x NCAR, por referencia                  *
 *                    cent   centroides, una matriz de tamano NGRUPOS x NCAR, por referencia              *
 *          Salida:   popul  grupo mas cercano a cada elemento, vector de tamano MAXE, por referencia      *
 **********************************************************************************************************/
void grupo_cercano (int nelem, float elem[][NCAR], float cent[][NCAR], int *popul) {
    int ngrupo;
    double adis, dmin;
    for (int i = 0; i < nelem; i++) {
        dmin = DBL_MAX;
        for (int j = 0; j < NGRUPOS; j++) {
            adis = gendist(elem[i], cent[j]); // elem[i] o &elem[i][0]
            if (adis < dmin) {
                dmin = adis;
                ngrupo = j;
            }
        }
        popul[i] = ngrupo;
    }
}
/**********************************************************************************************************
 * 3 - Funcion para calcular la densidad del grupo (dist. media entre todos sus elementos)                *
 *          Entrada:  elem     elementos, una matriz de tamano MAXE x NCAR, por referencia                *
 *                    listag   vector de NGRUPOS structs (informacion de grupos generados), por ref.      *
 *          Salida:   densidad densidad de los grupos (vector de tamano NGRUPOS, por referencia)          *
 **********************************************************************************************************/
void calcular_densidad (float elem[][NCAR], struct lista_grupos *listag, float *densidad) {
    for (int i = 0; i < NGRUPOS; i++) {
        int nelem = listag[i].nelemg;
        if (nelem < 2) {
            densidad[i] = 0;
        }
        else {
            int actg;
            double acum = 0.0, cont = 0.0;
            for (int j = 0; j < nelem; j++) {
                actg = listag[i].elemg[j];
                for (int k = j+1; k < nelem; k++) {
                    int othg = listag[i].elemg[k];
                    acum += gendist(elem[actg], elem[othg]);
                    cont += 1.0;
                }
            }
            densidad[i] = (float) (acum/cont);
        }
    }
}
/**********************************************************************************************************
 * 4 - Funcion para relizar el analisis de enfermedades                                                   *
 *          Entrada:  listag   vector de NGRUPOS structs (informacion de grupos generados), por ref.      *
 *                    enf      enfermedades, una matriz de tamano MAXE x TENF, por referencia             *
 *          Salida:   prob_enf vector de TENF structs (informacion del analisis realizado), por ref.      *
 **********************************************************************************************************/
void analizar_enfermedades (struct lista_grupos *listag, float enf[][TENF], struct analisis *prob_enf) {
    for (int i = 0; i < TENF; i++) {
        float mediamin = FLT_MAX, mediamax = FLT_MIN;
        int gmax, gmin;
        for (int j = 0; j < NGRUPOS; j++) {
            int nelem = listag[j].nelemg;
            float acum = 0;
            for (int k = 0; k < nelem; k++) {
                int actg = listag[j].elemg[k];
                acum += enf[actg][i];
            }
            float mediaact = acum/nelem;
            if (mediaact < mediamin) {
                mediamin = mediaact;
                gmin = j;
            }
            else if (mediaact >= mediamax) {
                mediamax = mediaact;
                gmax = j;
            }
        }
        prob_enf[i].max = mediamax;
        prob_enf[i].min = mediamin;
        prob_enf[i].gmax = gmax;
        prob_enf[i].gmin = gmin;
    }
}
