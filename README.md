# Genetics 🧬 💻
Programming project of the subject Computer Architecture given in second year of Computer Science Engineering in UPV/EHU.

### Warning!!! ⚠️
The execute may not be the same as the given ones in res1000.out and res.out. 
We do use a random number generator function and the fucntioning of this may depend on the compiler that you are using.

## Serial version (vserie)  
We developed the serial version of the program and some scripts to execute the program faster.

#### Code
The main functions that we have code are:
- **gendist(*elem 1, *elem2)**: Function to calculate vseries distance between two elements (Euclidean distance).
- **grupo_cercano(nelem, elem[][NCAR], cent[][NCAR], *popul)**: Function to calculate the closest cluster (closest centroid).
- **calcular_densidad(elem[][NCAR], *listag, *densidad)**: Function to calculate the density of the group (average distance between all its elements).
- **analizar_enfermedades (*listag, enf[][TENF], *prob_enf)**: Function to perform disease analysis.


## Parallel version (vparallel)
This version uses the API OpenMP. Is 20 times faster than the serial version, because we are distributing some tasks in different threads that can be computed simultaneously.
