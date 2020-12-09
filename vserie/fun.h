/**********************************************************************************************************
 *                      Cabeceras de las funciones utilizadas en el modulo gengrupos                      *
 *                                              fun.h                                                     *
 **********************************************************************************************************/

extern double gendist (float *elem1, float *elem2);
extern void grupo_cercano (int nelem, float elem[][NCAR], float cent[][NCAR], int *popul);
extern void calcular_densidad (float elem[][NCAR], struct lista_grupos *listag, float *densidad);
extern void analizar_enfermedades (struct lista_grupos *listag, float enf[][TENF], struct analisis *prob_enf);


