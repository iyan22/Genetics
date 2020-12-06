/**********************************************************************************************************
 *                              AC - OpenMP -- SERIE                                                      *
 *                  Rutinas que se utilizan en el modulo gengrupos_s.c                                    *
 *                                   fun_s.c                                                              *
 **********************************************************************************************************/
#include <math.h>
#include <float.h>

#include "defineg.h"    // Definiciones


/**********************************************************************************************************
 * 1 - Funcion para calcular la distancia vserie entre dos elementos (distancia euclidea)                 *
 *          Entrada:  2 elementos con NCAR caracteristicas (por referencia)                               *
 *          Salida:   distancia (double)                                                                  *
 **********************************************************************************************************/
double gendist (float *elem1, float *elem2) {
    float elem, acum = 0;
    for(int i = 0; i < NCAR; i++) {
        elem = *(elem1+(4*i)) - *(elem2+(4*i));
        elem = pow(elem, 2);
        acum += elem;
    }
    return sqrt(elem);
}
/**********************************************************************************************************
 * 2 - Funcion para calcular el grupo (cluster) mas cercano (centroide mas cercano)                       *
 *          Entrada:  nelem  numero de elementos, int                                                     *
 *                    elem   elementos, una matriz de tamano MAXE x NCAR, por referencia                  *
 *                    cent   centroides, una matriz de tamano NGRUPOS x NCAR, por referencia              *
 *          Salida:   popul  grupo mas cercano a cada elemento, vector de tamano MAXE, por referencia      *
 **********************************************************************************************************/
void grupo_cercano (int nelem, float elem[][NCAR], float cent[][NCAR], int *popul) {
    float aelem, adis, ngrupo, dmin = FLT_MAX;
    for (int i = 0; i < MAXE; i++) {
        for (int j = 0; j < NGRUPOS; j++) {
            adis = gendist((&elem)*i, (&cent)*j);
            if (adis < dmin) {
                dmin = adis;
                ngrupo = j;
            }
        }
        *(popul+(i*4)) = ngrupo;
    }
}
/**********************************************************************************************************
 * 3 - Funcion para calcular la densidad del grupo (dist. media entre todos sus elementos)                *
 *          Entrada:  elem     elementos, una matriz de tamano MAXE x NCAR, por referencia                *
 *                    listag   vector de NGRUPOS structs (informacion de grupos generados), por ref.      *
 *          Salida:   densidad densidad de los grupos (vector de tamano NGRUPOS, por referencia)          *
 **********************************************************************************************************/
void calcular_densidad (float elem[][NCAR], struct lista_grupos *listag, float *densidad) {
    int nelem = listag->nelemg;
    if (nelem < 2) {
        *densidad = 0;
    }
    else {
        float acum, actg;
        for (int i = 0; i < nelem; i++) {
            acum = 0;
            actg = listag->elemg[i];
            for (int j = 0; j < MAXE; j++) {
                if (j != actg) {
                    acum += gendist((&elem)*actg, (&elem)*j)
                }
            }
            *(densidad+(i*4)) = acum/MAXE-1;
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
    // TODO Realizar el analisis de enfermedades en los grupos:
    //        maximo y grupo en el que se da el maximo (para cada enfermedad)
    //        minimo y grupo en el que se da el minimo (para cada enfermedad)

}



