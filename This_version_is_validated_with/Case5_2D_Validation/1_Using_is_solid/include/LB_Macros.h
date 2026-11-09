
#ifndef LB_MACROS_H_
#define LB_MACROS_H_

//Logical constants
#define TRUE        1
#define FALSE       -1

//Macro functions
#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

// Delete 3D double array
inline void delete3D(double*** &arr, int dim1, int dim2) {// inline to avoid multiple definitions
    if (!arr) return;
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            delete[] arr[i][j];
        }
        delete[] arr[i];
    }
    delete[] arr;
    arr = nullptr;
}

// Delete 2D double array
inline void delete2D(double** &arr, int dim1) {
    if (!arr) return;
    for (int i = 0; i < dim1; ++i)
        delete[] arr[i];
    delete[] arr;
    arr = nullptr;
}

// Delete 2D int array
inline void delete2D_int(int** &arr, int dim1) {
    if (!arr) return;
    for (int i = 0; i < dim1; ++i)
        delete[] arr[i];
    delete[] arr;
    arr = nullptr;
}

#endif
