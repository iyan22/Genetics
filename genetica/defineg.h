/*************************************************************
  Definiciones utilizadas en los modulos de la aplicacion
                         defineg.h
**************************************************************/

#define MAXE     230000   // Numero de elementos (muestras)
#define NGRUPOS  100      // Numero de clusters
#define NCAR     40       // Dimensiones de cada muestra
#define TENF     20       // Tipos de enfermedad

#define DELTA    0.01     // Convergencia: cambio minimo en un centroide
#define MAXIT    10000    // Convergencia: numero de iteraciones maximo

/************************
 * Estructuras de datos *
 ************************/

// Informacion de los clusters
struct lista_grupos {
    int elemg[MAXE];   // Indices de los elementos
    int nelemg;        // Numero de elementos
};

// Resultados del analisis de enfermedades
struct analisis {
    float max, min;    // Maximo y minimo de todos los grupos
    int gmax, gmin;    // Grupos con los valores maximo y minimo
};