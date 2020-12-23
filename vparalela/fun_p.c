/**********************************************************************************************************
 *                              AC - OpenMP -- PARALELA                                                      *
 *                  Rutinas que se utilizan en el modulo gengrupos_s.c                                    *
 *                                   fun_s.c                                                              *
 **********************************************************************************************************/
#include <math.h>
#include <float.h>

#include "defineg.h"    // Definiciones


/**********************************************************************************************************
 * 1 - Funcion para calcular la distancia vserie entre dos elementos (distancia euclidea)                 *
 *          Entrada:  2 elementos con NCAR caracteristicas (por referencia)                               *
 *          Salida:  distancia (double)                                                                   *
 **********************************************************************************************************/
double gendist (float *elem1, float *elem2) {
    // TODO Calcular la distancia euclidea entre dos vectores
}

/**********************************************************************************************************
 * 2 - Funcion para calcular el grupo (cluster) mas cercano (centroide mas cercano)                       *
 *          Entrada:  nelem  numero de elementos, int                                                     *
 *                    elem   elementos, una matriz de tamano MAXE x NCAR, por referencia                  *
 *                    cent   centroides, una matriz de tamano NGRUPOS x NCAR, por referencia              *
 *          Salida:   popul  grupo mas ercano a cada elemento, vector de tamano MAXE, por referencia      *
 **********************************************************************************************************/
void grupo_cercano (int nelem, float elem[][NCAR], float cent[][NCAR], int *popul) {
    // TODO popul: grupo mas cercano a cada elemento
}

/**********************************************************************************************************
 * 3 - Funcion para calcular la densidad del grupo (dist. media entre todos sus elementos)                *
 *          Entrada:  elem     elementos, una matriz de tamano MAXE x NCAR, por referencia                *
 *                    listag   vector de NGRUPOS structs (informacion de grupos generados), por ref.      *
 *          Salida:   densidad densidad de los grupos (vector de tamano NGRUPOS, por referencia)          *
 **********************************************************************************************************/
void calcular_densidad (float elem[][NCAR], struct lista_grupos *listag, float *densidad) {
    // TODO Calcular la densidad de los grupos:
    //        media de las distancia entre todos los elementos del grupo
    //        si el numero de elementos del grupo es 0 o 1, densidad = 0
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


