/******************************************************************************/
/*                                                                            */
/**********			        TRIANGULACION DE DELAUNAY		         **********/
/*                                                                            */
/******************************************************************************/

/* Librerías                                                                  */
# define _USE_MATH_DEFINES
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <stdbool.h>
# include <math.h>
# include <float.h>
# include <stdbool.h>

/*********                    Estructuras de Datos                   **********/
/**                                                                          **/

/* Estructura para manejar las coordenadas de los puntos                      */
struct Punto {
    double x;
    double y;
    int indice;
};

struct Triangulo {
    struct Punto *vertices[3];
    struct Triangulo *vecinos[3];
    int indices[3];
    int esTrianguloSuper;
    int aristasRestringidas[3];
};

struct ListaTriangulos {
    struct Triangulo **triangulos;
    int numTriangulos;
    int capacidad;
};

struct Triangulacion {
    struct Triangulo *triangulos;  // Arreglo de triángulos
    int numTriangulos;             // Numero actual de triángulos
    int maxTriangulos;             // Capacidad máxima de triángulos
    struct Punto *puntos;          // Arreglo de puntos
    int numPuntos;                 // Numero de puntos
    int maxPuntos;                 // Capacidad máxima de puntos
    int numPuntosRegion1;          // Numero de puntos en la región 1
    // Función para agregar nuevos puntos
    void (*agregarPunto)(struct Triangulacion*, struct Punto*);
    struct Borde **bordes;
    int numBordes;
    int maxBordes;
};

struct Borde {
    struct Punto *p1;
    struct Punto *p2;
};

struct ColaRefinamiento {
    struct Triangulo **triangulos;
    struct Borde **bordes;
    int numTriangulos;
    int numBordes;
    int capacidad;
};

// Pool de memoria para gestión eficiente
struct PoolMemoria {
    void **elementos;
    int capacidadMaxima;
    int elementosLibres;
    int tamañoElemento;
    void *elementoActual;
};

// Estructura para manejar la entrada del archivo poly
struct EntradaPoly {
    // Vértices
    int numVertices;
    int numAtributos;
    int numMarcadores;
    struct Punto *vertices;
    double *atributos;  // Atributos por vértice
    int *marcadores;    // Marcadores de frontera

    // Segmentos
    int numSegmentos;
    struct Segmento *segmentos;
    int *marcadoresSegmentos;

    // Agujeros
    int numAgujeros;
    struct Punto *agujeros;

    // Regiones
    int numRegiones;
    struct Region *regiones;
};

struct Region {
    double x, y;          // Punto dentro de la región
    int atributo;         // Atributo de la región
    double areaMaxima;    // Restricción de área
};

struct Segmento {
    int v1, v2;           // Índices de los vértices que forman el segmento
    int marcador;         // Marcador del segmento (ej: frontera)
};


/*********                    Prototipos de Funciones                **********/
/**                                                                          **/
static int leerLineaNoVacia(FILE* archivo, char* linea, int maxLen);
static int leerSegmentos(FILE* archivo, struct EntradaPoly* entrada);
static int leerAgujeros(FILE* archivo, struct EntradaPoly* entrada);
static int leerRegiones(FILE* archivo, struct EntradaPoly* entrada);
struct EntradaPoly* leerArchivoPoly(const char *nombreArchivo);
void guardarArchivoNode(struct Triangulacion* tr, const char* nombreArchivo);
void guardarArchivoEle(struct Triangulacion *tr, const char *nombreArchivo);
struct PoolMemoria* inicializarPoolMemoria(int capacidadMaxima, int tamañoElemento);
void* reservarMemoria(struct PoolMemoria* pool);
void liberarMemoria(struct PoolMemoria* pool, void* elemento);
void liberarEntradaPoly(struct EntradaPoly* entrada);
double areaTriangulo(struct Punto *a, struct Punto *b, struct Punto *c);
void intercambiarDiagonal(struct Triangulacion *tr, struct Triangulo *t1, struct Triangulo *t2);
double determinante3x3(double matriz[3][3]);
struct Punto* calcularCircuncentro(struct Triangulo *t);
double calcularAreaTriangulo(struct Triangulo *t);
double calcularAngulo(struct Punto *p1, struct Punto *p2, struct Punto *p3);
double distanciaEntrePuntos(struct Punto *p1, struct Punto *p2);
int compararPuntos(const void *a, const void *b);
bool puntoEnTriangulo(struct Triangulo *t, struct Punto *p);
int tieneVerticeArtificial(struct Triangulo *t);
int tieneVertice(struct Triangulo *t, struct Punto *p);
int puntoEnCircunscrito(struct Triangulo *t, struct Punto *p);
bool necesitaRefinamiento(struct Triangulo *t, double anguloMinimo, double areaMaxima);
bool esPuntoCercaDelBorde(struct Punto *p, struct Triangulacion *tr);
struct Punto* calcularPuntoInterseccion(struct Triangulo *t1, struct Triangulo *t2,
                                       struct Punto *p1, struct Punto *p2);
bool encontrarInterseccion(double x1, double y1, double x2, double y2,
                          double x3, double y3, double x4, double y4,
                          double *x, double *y);
bool esPuntoExtremo(double x, double y, struct Punto *p1, struct Punto *p2);
int encontrarBorde(struct Triangulo *t, struct Punto *p1, struct Punto *p2);
bool necesitaIntercambio(struct Triangulo *t1, struct Triangulo *t2);
void actualizarVecinos(struct Triangulacion *tr);
void agregarTrianguloACola(struct ColaRefinamiento *cola, struct Triangulo *t);
struct Triangulo* extraerTriangulo(struct ColaRefinamiento *cola);
bool dentroLimites(struct Punto *p, struct Triangulacion *tr);
struct Triangulacion* inicializarTriangulacion(struct Punto* puntos, int numPuntos, int numPuntosRegion1);
struct Triangulo* crearTriangulo(struct Punto *v1, struct Punto *v2, struct Punto *v3);
struct Triangulo* crearSuperTriangulo(struct Punto *puntos, int numPuntos);
void agregarTriangulo(struct Triangulacion *tr, struct Triangulo *t);
void eliminarTriangulo(struct Triangulacion *tr, struct Triangulo *t);
void dividirTriangulo(struct Triangulacion *tr, struct Triangulo *t, struct Punto *p);
void eliminarTriangulosSuper(struct Triangulacion *tr);
void liberarTriangulacion(struct Triangulacion *tr);
struct ListaTriangulos* encontrarTriangulosIntersectados(struct Triangulacion *tr, 
                                                        struct Punto *p1, 
                                                        struct Punto *p2);
void divideVencerasDelaunay(struct Triangulacion *tr, int inicio, int fin);
void combinarTriangulaciones(struct Triangulacion *tr, int inicio, int medio, int fin);
struct Borde* encontrarBaseInicial(struct Triangulacion *tr, int medio);
struct Punto* encontrarCandidatoIzquierda(struct Triangulacion *tr, struct Borde *base);
struct Punto* encontrarCandidatoDerecha(struct Triangulacion *tr, struct Borde *base);
bool puntoEnCircunferencia(struct Punto *p1, struct Punto *p2, struct Punto *p3, struct Punto *punto);
double calcularAngulo(struct Punto *p1, struct Punto *p2, struct Punto *p3);
int compararPuntosX(const void *a, const void *b);
void crearSegmento(struct Triangulacion *tr, struct Punto *p1, struct Punto *p2);
struct Triangulacion* triangulacionDelaunay(struct Punto* puntos, int numPuntos, int numPuntosRegion1);
struct ColaRefinamiento* inicializarColaRefinamiento(int capacidad);
void agregarPuntoATriangulacion(struct Triangulacion *tr, struct Punto *p);
void agregarTrianguloACola(struct ColaRefinamiento *cola, struct Triangulo *t);
struct Triangulo* extraerTriangulo(struct ColaRefinamiento *cola);
bool dentroLimites(struct Punto *p, struct Triangulacion *tr);
void refinarMalla(struct Triangulacion *tr, double anguloMinimo, double areaMaxima);
void liberarColaRefinamiento(struct ColaRefinamiento *cola);
struct Punto* calcularCircuncentro(struct Triangulo *t);
bool estaDentroDeLimites(struct Triangulacion *tr, struct Punto *p);
bool hayPuntoCercano(struct Triangulacion *tr, struct Punto *p);
double calcularAreaMaximaPermitida(struct Triangulacion *tr);
bool necesitaRefinamiento(struct Triangulo *t, double anguloMinimo, double areaMaxima);
void calcularAngulos(struct Triangulo *t, double *angulos);
double calcularAreaTriangulo2(struct Punto *p1, struct Punto *p2, struct Punto *p3);
struct Borde* crearBorde(struct Triangulacion *tr, struct Punto *p1, struct Punto *p2);
void agregarBorde(struct Triangulacion *tr, struct Punto *p1, struct Punto *p2);
void combinarTriangulacionesOptimizado(struct Triangulacion *tr, int inicio, int medio, int fin);
bool existeTriangulo(struct Triangulacion *tr, struct Punto *p1, struct Punto *p2, struct Punto *p3);
void menu(void);
void info(void);

/*********                         Funciones                         **********/
/**                                                                          **/


/* Funciones de lectura y de gestión de archivos                              */

static int leerLineaNoVacia(FILE* archivo, char* linea, int maxLen) {
    while (fgets(linea, maxLen, archivo)) {
        // Saltar comentarios y líneas vacías
        if (linea[0] != '#' && linea[0] != '\n' && linea[0] != '\0') {
            return 1;
        }
    }
    return 0;
}

static int leerSegmentos(FILE* archivo, struct EntradaPoly* entrada) {
    char linea[256];
    
    // Leer cabecera de segmentos
    if (!leerLineaNoVacia(archivo, linea, sizeof(linea))) {
        return 0;
    }

    int numMarcadores;
    if (sscanf(linea, "%d %d", &entrada->numSegmentos, &numMarcadores) != 2) {
        return 0;
    }

    if (entrada->numSegmentos > 0) {
        entrada->segmentos = malloc(entrada->numSegmentos * sizeof(struct Segmento));
        if (numMarcadores) {
            entrada->marcadoresSegmentos = malloc(entrada->numSegmentos * sizeof(int));
        }

        for (int i = 0; i < entrada->numSegmentos; i++) {
            if (!leerLineaNoVacia(archivo, linea, sizeof(linea))) {
                return 0;
            }

            int indice;
            if (numMarcadores) {
                if (sscanf(linea, "%d %d %d %d", &indice, 
                          &entrada->segmentos[i].v1, 
                          &entrada->segmentos[i].v2,
                          &entrada->segmentos[i].marcador) != 4) {
                    return 0;
                }
            } else {
                if (sscanf(linea, "%d %d %d", &indice, 
                          &entrada->segmentos[i].v1, 
                          &entrada->segmentos[i].v2) != 3) {
                    return 0;
                }
                entrada->segmentos[i].marcador = 0;
            }
        }
    }
    return 1;
}

static int leerAgujeros(FILE* archivo, struct EntradaPoly* entrada) {
    char linea[256];
    
    if (!leerLineaNoVacia(archivo, linea, sizeof(linea))) {
        return 0;
    }

    if (sscanf(linea, "%d", &entrada->numAgujeros) != 1) {
        return 0;
    }

    if (entrada->numAgujeros > 0) {
        entrada->agujeros = malloc(entrada->numAgujeros * sizeof(struct Punto));
        
        for (int i = 0; i < entrada->numAgujeros; i++) {
            if (!leerLineaNoVacia(archivo, linea, sizeof(linea))) {
                return 0;
            }

            int indice;
            if (sscanf(linea, "%d %lf %lf", &indice, 
                      &entrada->agujeros[i].x, 
                      &entrada->agujeros[i].y) != 3) {
                return 0;
            }
            entrada->agujeros[i].indice = indice;
        }
    }
    return 1;
}

static int leerRegiones(FILE* archivo, struct EntradaPoly* entrada) {
    char linea[256];
    
    if (!leerLineaNoVacia(archivo, linea, sizeof(linea))) {
        return 0;
    }

    if (sscanf(linea, "%d", &entrada->numRegiones) != 1) {
        return 0;
    }

    if (entrada->numRegiones > 0) {
        entrada->regiones = malloc(entrada->numRegiones * sizeof(struct Region));
        
        for (int i = 0; i < entrada->numRegiones; i++) {
            if (!leerLineaNoVacia(archivo, linea, sizeof(linea))) {
                return 0;
            }

            int indice;
            if (sscanf(linea, "%d %lf %lf %d %lf", &indice, 
                      &entrada->regiones[i].x, 
                      &entrada->regiones[i].y,
                      &entrada->regiones[i].atributo,
                      &entrada->regiones[i].areaMaxima) != 5) {
                return 0;
            }
        }
    }
    return 1;
}

struct EntradaPoly* leerArchivoPoly(const char *nombreArchivo) {
    FILE *archivo = fopen(nombreArchivo, "r");
    if (!archivo) {
        printf("[ERROR] No se pudo abrir el archivo %s\n", nombreArchivo);
        return NULL;
    }

    struct EntradaPoly *entrada = (struct EntradaPoly*)malloc(sizeof(struct EntradaPoly));
    if (!entrada) {
        printf("[ERROR] No se pudo asignar memoria para la estructura de entrada\n");
        fclose(archivo);
        return NULL;
    }

    // Inicializar con valores por defecto
    memset(entrada, 0, sizeof(struct EntradaPoly));

    char linea[256];
    int dim; // Añadimos la variable dim aquí

    // Leer sección de vértices
    if (!leerLineaNoVacia(archivo, linea, sizeof(linea))) {
        printf("[ERROR] Archivo vacío o mal formateado\n");
        liberarEntradaPoly(entrada);
        fclose(archivo);
        return NULL;
    }

    // Procesar primera línea: número de vértices, dimensión, atributos, marcadores
    if (sscanf(linea, "%d %d %d %d", &entrada->numVertices, &dim, 
               &entrada->numAtributos, &entrada->numMarcadores) != 4) {
        printf("[ERROR] Error al leer la cabecera de vértices\n");
        liberarEntradaPoly(entrada);
        fclose(archivo);
        return NULL;
    }

    // Verificar dimensión
    if (dim != 2) {
        printf("[ERROR] Solo se soportan archivos 2D (dimensión = %d)\n", dim);
        liberarEntradaPoly(entrada);
        fclose(archivo);
        return NULL;
    }

    // Asignar memoria para vértices y sus atributos
    entrada->vertices = malloc(entrada->numVertices * sizeof(struct Punto));
    if (entrada->numAtributos > 0) {
        entrada->atributos = malloc(entrada->numVertices * entrada->numAtributos * sizeof(double));
    }
    if (entrada->numMarcadores > 0) {
        entrada->marcadores = malloc(entrada->numVertices * sizeof(int));
    }

    // Leer los vértices
    for (int i = 0; i < entrada->numVertices; i++) {
        if (!leerLineaNoVacia(archivo, linea, sizeof(linea))) {
            printf("[ERROR] Error al leer vértice %d\n", i);
            liberarEntradaPoly(entrada);
            fclose(archivo);
            return NULL;
        }

        int indice;
        double x, y;
        if (sscanf(linea, "%d %lf %lf", &indice, &x, &y) != 3) {
            printf("[ERROR] Formato inválido para vértice %d\n", i);
            liberarEntradaPoly(entrada);
            fclose(archivo);
            return NULL;
        }

        entrada->vertices[i].x = x;
        entrada->vertices[i].y = y;
        entrada->vertices[i].indice = indice - 1; // Convertir a base-0

        // Leer atributos si existen
        char *token = strtok(linea, " \t");
        for (int j = 0; j < 3; j++) token = strtok(NULL, " \t"); // Saltar índice y coordenadas
        for (int j = 0; j < entrada->numAtributos; j++) {
            token = strtok(NULL, " \t");
            if (token) {
                entrada->atributos[i * entrada->numAtributos + j] = atof(token);
            }
        }

        // Leer marcador si existe
        if (entrada->numMarcadores > 0) {
            token = strtok(NULL, " \t");
            if (token) {
                entrada->marcadores[i] = atoi(token);
            }
        }
    }

    // Leer segmentos
    if (!leerSegmentos(archivo, entrada)) {
        printf("[ERROR] Error al leer segmentos\n");
        liberarEntradaPoly(entrada);
        fclose(archivo);
        return NULL;
    }

    // Leer agujeros
    if (!leerAgujeros(archivo, entrada)) {
        printf("[ERROR] Error al leer agujeros\n");
        liberarEntradaPoly(entrada);
        fclose(archivo);
        return NULL;
    }

    // Leer regiones
    if (!leerRegiones(archivo, entrada)) {
        printf("[ERROR] Error al leer regiones\n");
        liberarEntradaPoly(entrada);
        fclose(archivo);
        return NULL;
    }

    fclose(archivo);
    return entrada;
}

void guardarArchivoNode(struct Triangulacion* tr, const char* nombreArchivo) {
    FILE* archivo = fopen(nombreArchivo, "w");
    if (!archivo) {
        printf("[ERROR] No se pudo crear el archivo %s\n", nombreArchivo);
        return;
    }

    // Escribir encabezado
    fprintf(archivo, "%d 2 0 1\n", tr->numPuntos);
    
    // Guardar los puntos con sus coordenadas originales
    for(int i = 0; i < tr->numPuntos; i++) {
        int region;
        // Los primeros n1 puntos son de la región 1, el resto de la región 2
        if (i < tr->numPuntosRegion1) {
            region = 1;
        } else {
            region = 2;
        }
        
        fprintf(archivo, "%d %.6f %.6f %d\n", 
                i + 1,                // índice del punto
                tr->puntos[i].x,     // coordenada x
                tr->puntos[i].y,     // coordenada y
                region);             // región (1 o 2)
    }

    fprintf(archivo, "# Generated by Delaunay Triangulation\n");
    fclose(archivo);
}

void guardarArchivoEle(struct Triangulacion *tr, const char *nombreArchivo) {
    printf("Iniciando guardarArchivoEle...\n");
    
    if (tr == NULL) {
        printf("Error: tr es NULL\n");
        return;
    }

    printf("Número de triángulos en tr: %d\n", tr->numTriangulos);
    
    FILE *archivo = fopen(nombreArchivo, "w");
    if (archivo == NULL) {
        printf("Error: No se puede crear el archivo %s. errno: %d\n", nombreArchivo, errno);
        return;
    }

    // Contar solo triángulos válidos y verificar índices
    int triangulos_validos = 0;
    int i;
    for (i = 0; i < tr->numTriangulos; i++) {
        if (!tr->triangulos[i].esTrianguloSuper) {
            // Verificar que los índices sean válidos
            if (tr->triangulos[i].vertices[0]->indice >= 0 && 
                tr->triangulos[i].vertices[0]->indice < tr->numPuntos &&
                tr->triangulos[i].vertices[1]->indice >= 0 && 
                tr->triangulos[i].vertices[1]->indice < tr->numPuntos &&
                tr->triangulos[i].vertices[2]->indice >= 0 && 
                tr->triangulos[i].vertices[2]->indice < tr->numPuntos) {
                triangulos_validos++;
            } else {
                printf("Triángulo %d ignorado por índices inválidos: %d %d %d\n",
                    i,
                    tr->triangulos[i].vertices[0]->indice,
                    tr->triangulos[i].vertices[1]->indice,
                    tr->triangulos[i].vertices[2]->indice);
            }
        }
    }
    
    printf("Número de triángulos válidos: %d\n", triangulos_validos);
    fprintf(archivo, "%d  3  0\n", triangulos_validos);
    
    int indice = 1;
    for (i = 0; i < tr->numTriangulos; i++) {
        if (!tr->triangulos[i].esTrianguloSuper) {
            // Solo escribir triángulos con índices válidos
            if (tr->triangulos[i].vertices[0]->indice >= 0 && 
                tr->triangulos[i].vertices[0]->indice < tr->numPuntos &&
                tr->triangulos[i].vertices[1]->indice >= 0 && 
                tr->triangulos[i].vertices[1]->indice < tr->numPuntos &&
                tr->triangulos[i].vertices[2]->indice >= 0 && 
                tr->triangulos[i].vertices[2]->indice < tr->numPuntos) {
                
                fprintf(archivo, "%d  %d  %d  %d\n",
                        indice++,
                        tr->triangulos[i].vertices[0]->indice + 1,  // Convertir a base-1 para MATLAB
                        tr->triangulos[i].vertices[1]->indice + 1,
                        tr->triangulos[i].vertices[2]->indice + 1);
                
                printf("Escribiendo triángulo %d: (%d, %d, %d)\n", 
                       indice - 1,
                       tr->triangulos[i].vertices[0]->indice + 1,
                       tr->triangulos[i].vertices[1]->indice + 1,
                       tr->triangulos[i].vertices[2]->indice + 1);
            }
        }
    }

    fprintf(archivo, "# Generated by Delaunay Triangulation\n");
    fclose(archivo);
    printf("Archivo .ele guardado exitosamente.\n");
}


/* Funciones de Gestión de Memoria                                                */

struct PoolMemoria* inicializarPool(int capacidadMaxima, int tamañoElemento) {
    struct PoolMemoria* pool = (struct PoolMemoria*)malloc(sizeof(struct PoolMemoria));
    if (!pool) {
        printf("[ERROR] No se pudo crear el pool de memoria\n");
        return NULL;
    }

    pool->elementos = (void**)malloc(capacidadMaxima * sizeof(void*));
    if (!pool->elementos) {
        printf("[ERROR] No se pudo asignar memoria para los elementos del pool\n");
        free(pool);
        return NULL;
    }

    pool->elementoActual = malloc(capacidadMaxima * tamañoElemento);
    if (!pool->elementoActual) {
        printf("[ERROR] No se pudo asignar el bloque de memoria principal\n");
        free(pool->elementos);
        free(pool);
        return NULL;
    }

    pool->capacidadMaxima = capacidadMaxima;
    pool->elementosLibres = capacidadMaxima;
    pool->tamañoElemento = tamañoElemento;

    // Inicializar la lista de elementos libres
    for (int i = 0; i < capacidadMaxima; i++) {
        pool->elementos[i] = (char*)pool->elementoActual + (i * tamañoElemento);
    }

    return pool;
}

void* obtenerDelPool(struct PoolMemoria* pool) {
    if (!pool || pool->elementosLibres <= 0) {
        return NULL;
    }

    void* elemento = pool->elementos[--pool->elementosLibres];
    memset(elemento, 0, pool->tamañoElemento);
    return elemento;
}

void devolverAlPool(struct PoolMemoria* pool, void* elemento) {
    if (!pool || !elemento || pool->elementosLibres >= pool->capacidadMaxima) {
        return;
    }

    pool->elementos[pool->elementosLibres++] = elemento;
}

void liberarPool(struct PoolMemoria* pool) {
    if (pool) {
        free(pool->elementoActual);
        free(pool->elementos);
        free(pool);
    }
}

void liberarEntradaPoly(struct EntradaPoly* entrada) {
    if (entrada) {
        free(entrada->vertices);
        free(entrada->atributos);
        free(entrada->marcadores);
        free(entrada->segmentos);
        free(entrada->marcadoresSegmentos);
        free(entrada->agujeros);
        free(entrada->regiones);
        free(entrada);
    }
}


/* Funciones de Operaciones Aritméticas                                           */

// Función para calcular el área de un triángulo usando determinantes
double areaTriangulo(struct Punto *a, struct Punto *b, struct Punto *c) {
    return fabs((b->x - a->x) * (c->y - a->y) - (c->x - a->x) * (b->y - a->y)) / 2.0;
}

void intercambiarDiagonal(struct Triangulacion *tr, struct Triangulo *t1, struct Triangulo *t2) {
    // Implementación básica de intercambio de diagonal
    struct Punto *p1 = NULL, *p2 = NULL, *p3 = NULL, *p4 = NULL;
    
    // Encontrar los cuatro puntos del cuadrilátero
    for(int i = 0; i < 3; i++) {
        if(!tieneVertice(t2, t1->vertices[i])) p1 = t1->vertices[i];
        if(!tieneVertice(t1, t2->vertices[i])) p3 = t2->vertices[i];
    }
    
    // Encontrar los puntos compartidos
    for(int i = 0; i < 3; i++) {
        if(t1->vertices[i] != p1) {
            if(p2 == NULL) p2 = t1->vertices[i];
            else p4 = t1->vertices[i];
        }
    }
    
    // Crear nuevos triángulos
    *t1 = *crearTriangulo(p1, p3, p2);
    *t2 = *crearTriangulo(p1, p3, p4);
}

// Función para calcular el determinante de una matriz 3x3
double determinante3x3(double matriz[3][3]) {
    return matriz[0][0] * (matriz[1][1] * matriz[2][2] - matriz[1][2] * matriz[2][1])
         - matriz[0][1] * (matriz[1][0] * matriz[2][2] - matriz[1][2] * matriz[2][0])
         + matriz[0][2] * (matriz[1][0] * matriz[2][1] - matriz[1][1] * matriz[2][0]);
}

struct Punto* calcularCircuncentro(struct Triangulo *t) {
    struct Punto *p1 = t->vertices[0];
    struct Punto *p2 = t->vertices[1];
    struct Punto *p3 = t->vertices[2];
    
    double d = 2 * (p1->x * (p2->y - p3->y) + 
                    p2->x * (p3->y - p1->y) + 
                    p3->x * (p1->y - p2->y));
    
    if (fabs(d) < 1e-10) return NULL;  // Triángulo degenerado
    
    struct Punto *c = malloc(sizeof(struct Punto));
    if (!c) return NULL;
    
    c->x = ((p1->x*p1->x + p1->y*p1->y) * (p2->y - p3->y) +
            (p2->x*p2->x + p2->y*p2->y) * (p3->y - p1->y) +
            (p3->x*p3->x + p3->y*p3->y) * (p1->y - p2->y)) / d;
            
    c->y = ((p1->x*p1->x + p1->y*p1->y) * (p3->x - p2->x) +
            (p2->x*p2->x + p2->y*p2->y) * (p1->x - p3->x) +
            (p3->x*p3->x + p3->y*p3->y) * (p2->x - p1->x)) / d;
            
    return c;
}

double calcularAreaTriangulo(struct Triangulo *t) {
    if (t == NULL) return 0.0;

    struct Punto *p1 = t->vertices[0];
    struct Punto *p2 = t->vertices[1];
    struct Punto *p3 = t->vertices[2];

    // Fórmula del área usando determinante
    double area = fabs((p1->x * (p2->y - p3->y) + 
                       p2->x * (p3->y - p1->y) + 
                       p3->x * (p1->y - p2->y)) / 2.0);

    return area;
}

double calcularAngulo(struct Punto *p1, struct Punto *p2, struct Punto *p3) {
    double dx1 = p2->x - p1->x;
    double dy1 = p2->y - p1->y;
    double dx2 = p3->x - p1->x;
    double dy2 = p3->y - p1->y;
    
    // Calcular el ángulo usando el producto punto
    double dot = dx1*dx2 + dy1*dy2;
    double len1 = sqrt(dx1*dx1 + dy1*dy1);
    double len2 = sqrt(dx2*dx2 + dy2*dy2);
    
    if(len1 == 0 || len2 == 0) return -1;
    
    double cosAngulo = dot/(len1*len2);
    if(cosAngulo > 1) cosAngulo = 1;
    if(cosAngulo < -1) cosAngulo = -1;
    
    return acos(cosAngulo);
}

double distanciaEntrePuntos(struct Punto *p1, struct Punto *p2) {
    return sqrt(pow(p2->x - p1->x, 2) + pow(p2->y - p1->y, 2));
}


/* Funciones de verificación y comparación                                        */

// Función para comparar puntos (primero por x, luego por y)
int compararPuntos(const void *a, const void *b) {
    struct Punto *p1 = (struct Punto *)a;
    struct Punto *p2 = (struct Punto *)b;
    
    if (p1->x != p2->x) {
        return (p1->x > p2->x) ? 1 : -1;
    }
    return (p1->y > p2->y) ? 1 : -1;
}

// Función para verificar si un punto está dentro de un triángulo usando áreas
bool puntoEnTriangulo(struct Triangulo *t, struct Punto *p) {
    // Cálculo de áreas usando determinantes
    double area = areaTriangulo(t->vertices[0], t->vertices[1], t->vertices[2]);
    double area1 = areaTriangulo(p, t->vertices[1], t->vertices[2]);
    double area2 = areaTriangulo(t->vertices[0], p, t->vertices[2]);
    double area3 = areaTriangulo(t->vertices[0], t->vertices[1], p);

    // Usar tolerancia para comparaciones de punto flotante
    const double EPSILON = 1e-10;
    double total = area1 + area2 + area3;
    
    return fabs(area - total) < EPSILON;
}

// Función para verificar si un triángulo tiene vértices artificiales
int tieneVerticeArtificial(struct Triangulo *t) {
    return t->vertices[0]->indice < 0 || 
           t->vertices[1]->indice < 0 || 
           t->vertices[2]->indice < 0;
}

// Función auxiliar para verificar si un triángulo tiene un vértice
int tieneVertice(struct Triangulo *t, struct Punto *p) {
    if(t == NULL || p == NULL) return 0;
    
    return (t->vertices[0] == p || 
            t->vertices[1] == p || 
            t->vertices[2] == p);
}

// Función para verificar si un punto está dentro del círculo circunscrito
int puntoEnCircunscrito(struct Triangulo *t, struct Punto *p) {
    struct Punto *a = t->vertices[0];
    struct Punto *b = t->vertices[1];
    struct Punto *c = t->vertices[2];
    
    // Matriz para el test del círculo circunscrito
    double matriz[3][3] = {
        {a->x - p->x, a->y - p->y, (a->x - p->x) * (a->x - p->x) + (a->y - p->y) * (a->y - p->y)},
        {b->x - p->x, b->y - p->y, (b->x - p->x) * (b->x - p->x) + (b->y - p->y) * (b->y - p->y)},
        {c->x - p->x, c->y - p->y, (c->x - p->x) * (c->x - p->x) + (c->y - p->y) * (c->y - p->y)}
    };
    
    // El signo del determinante indica si el punto está dentro o fuera
    double det = determinante3x3(matriz);
    
    // Usamos una pequeña tolerancia para manejar errores de punto flotante
    double tolerancia = 1e-10;
    
    // Si el determinante es positivo, el punto está dentro del círculo
    return det > tolerancia;
}

bool necesitaRefinamiento(struct Triangulo *t, double anguloMinimo, double areaMaxima) {
    // Verificar área
    double area = calcularAreaTriangulo(t);
    if (area > areaMaxima) return true;
    
    // Verificar ángulos
    double angulos[3];
    calcularAngulos(t, angulos);
    
    for (int i = 0; i < 3; i++) {
        if (angulos[i] < anguloMinimo) return true;
    }
    
    return false;
}

void calcularAngulos(struct Triangulo *t, double *angulos) {
    for (int i = 0; i < 3; i++) {
        struct Punto *p1 = t->vertices[i];
        struct Punto *p2 = t->vertices[(i+1)%3];
        struct Punto *p3 = t->vertices[(i+2)%3];
        
        double dx1 = p2->x - p1->x;
        double dy1 = p2->y - p1->y;
        double dx2 = p3->x - p1->x;
        double dy2 = p3->y - p1->y;
        
        angulos[i] = fabs(atan2(dx1*dy2 - dy1*dx2, dx1*dx2 + dy1*dy2));
    }
}

bool esPuntoCercaDelBorde(struct Punto *p, struct Triangulacion *tr) {
    const double DIST_MIN_BORDE = 0.05;
    const int NUM_VECINOS = 2; // Reducimos a 2 para ser más selectivos
    int puntos_cercanos = 0;
    
    // Primero verificar si está muy lejos verticalmente
    double ymin = tr->puntos[0].y;
    double ymax = tr->puntos[0].y;
    for(int i = 1; i < tr->numPuntos; i++) {
        if(tr->puntos[i].y < ymin) ymin = tr->puntos[i].y;
        if(tr->puntos[i].y > ymax) ymax = tr->puntos[i].y;
    }
    
    if(p->y < ymin - 0.1 || p->y > ymax + 0.1) {
        return true;  // Rechazar puntos muy alejados verticalmente
    }
    
    // Verificar puntos cercanos
    for(int i = 0; i < tr->numPuntos; i++) {
        if(distanciaEntrePuntos(&tr->puntos[i], p) < DIST_MIN_BORDE) {
            puntos_cercanos++;
            if(puntos_cercanos >= NUM_VECINOS) {
                return true;
            }
        }
    }
    return false;
}

struct Punto* calcularPuntoInterseccion(struct Triangulo *t1, struct Triangulo *t2,
                                       struct Punto *p1, struct Punto *p2) {
    // Crear punto para almacenar la intersección
    struct Punto *interseccion = malloc(sizeof(struct Punto));
    if (!interseccion) return NULL;

    // Obtener los segmentos del triángulo que podrían intersectar
    for (int i = 0; i < 3; i++) {
        struct Punto *a = t1->vertices[i];
        struct Punto *b = t1->vertices[(i+1)%3];
        
        // Verificar si este segmento intersecta con el segmento restricción
        double x, y;
        if (encontrarInterseccion(p1->x, p1->y, p2->x, p2->y,
                                 a->x, a->y, b->x, b->y,
                                 &x, &y)) {
            // Verificar que el punto no sea uno de los extremos
            if (!esPuntoExtremo(x, y, p1, p2) && 
                !esPuntoExtremo(x, y, a, b)) {
                interseccion->x = x;
                interseccion->y = y;
                interseccion->indice = -1; // Punto temporal
                return interseccion;
            }
        }
    }
    
    free(interseccion);
    return NULL;
}

// Función auxiliar para encontrar la intersección entre dos segmentos
bool encontrarInterseccion(double x1, double y1, double x2, double y2,
                          double x3, double y3, double x4, double y4,
                          double *x, double *y) {
    // Calcular denominador
    double denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1);
    if (fabs(denom) < 1e-10) return false; // Segmentos paralelos

    // Calcular parámetros de intersección
    double ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom;
    double ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom;

    // Verificar si la intersección está dentro de ambos segmentos
    if (ua >= 0 && ua <= 1 && ub >= 0 && ub <= 1) {
        *x = x1 + ua*(x2-x1);
        *y = y1 + ua*(y2-y1);
        return true;
    }

    return false;
}

// Función auxiliar para verificar si un punto es extremo de un segmento
bool esPuntoExtremo(double x, double y, struct Punto *p1, struct Punto *p2) {
    const double EPSILON = 1e-10;
    return (fabs(x - p1->x) < EPSILON && fabs(y - p1->y) < EPSILON) ||
           (fabs(x - p2->x) < EPSILON && fabs(y - p2->y) < EPSILON);
}

bool esDelaunay(struct Punto *candidatoIzq, struct Punto *candidatoDer, struct Borde *base) {
    if (!candidatoIzq) return true;   // Si no hay candidato izquierdo, usar derecho
    if (!candidatoDer) return false;  // Si no hay candidato derecho, usar izquierdo
    
    // Verificar el criterio del círculo vacío
    return !puntoEnCircunferencia(candidatoIzq, base->p1, base->p2, candidatoDer);
}

double orientacion(struct Punto *p1, struct Punto *p2, struct Punto *p3) {
    return (p2->x - p1->x) * (p3->y - p1->y) - 
           (p3->x - p1->x) * (p2->y - p1->y);
}

bool puntoEnCircunferencia(struct Punto *p1, struct Punto *p2, 
                          struct Punto *p3, struct Punto *punto) {
    double ax = p2->x - p1->x;
    double ay = p2->y - p1->y;
    double bx = p3->x - p1->x;
    double by = p3->y - p1->y;
    
    double c = ax * ax + ay * ay;
    double d = bx * bx + by * by;
    double e = ax * by - ay * bx;
    
    if (e == 0) return false;  // Puntos colineales
    
    double px = punto->x - p1->x;
    double py = punto->y - p1->y;
    
    double det = (px * px + py * py) * (ax * by - ay * bx) -
                (c * (px * by - py * bx) + d * (py * ax - px * ay));
                
    return det > 0;
}


/* Funciones auxiliares                                                        */

bool necesitaIntercambio(struct Triangulo *t1, struct Triangulo *t2) {
    // Implementación del criterio de Delaunay
    struct Punto *p4 = NULL;
    for(int i = 0; i < 3; i++) {
        if(!tieneVertice(t1, t2->vertices[i])) {
            p4 = t2->vertices[i];
            break;
        }
    }
    return puntoEnCircunscrito(t1, p4);
}

// Función para encontrar el índice del vértice opuesto a un borde en un triángulo
int encontrarBorde(struct Triangulo *t, struct Punto *p1, struct Punto *p2) {
    int i;
    for(i = 0; i < 3; i++) {
        if((t->vertices[i] == p1 && t->vertices[(i+1)%3] == p2) ||
           (t->vertices[i] == p2 && t->vertices[(i+1)%3] == p1)) {
            return i;
        }
    }
    return -1;  // No se encontró el borde
}

// Función para actualizar las relaciones de vecindad
void actualizarVecinos(struct Triangulacion *tr) {
    int i, j, k;
    
    // Primero, limpiamos todas las relaciones de vecindad
    for(i = 0; i < tr->numTriangulos; i++) {
        tr->triangulos[i].vecinos[0] = NULL;
        tr->triangulos[i].vecinos[1] = NULL;
        tr->triangulos[i].vecinos[2] = NULL;
    }
    
    // Luego, buscamos los vecinos para cada triángulo
    for(i = 0; i < tr->numTriangulos; i++) {
        for(j = i + 1; j < tr->numTriangulos; j++) {
            // Para cada par de triángulos, verificamos si comparten un borde
            for(k = 0; k < 3; k++) {
                struct Punto *p1 = tr->triangulos[i].vertices[k];
                struct Punto *p2 = tr->triangulos[i].vertices[(k+1)%3];
                
                // Buscar este borde en el otro triángulo
                int bordeComun = encontrarBorde(&tr->triangulos[j], p1, p2);
                if(bordeComun != -1) {
                    // Los triángulos son vecinos
                    tr->triangulos[i].vecinos[k] = &tr->triangulos[j];
                    tr->triangulos[j].vecinos[bordeComun] = &tr->triangulos[i];
                }
            }
        }
    }
}

void agregarTrianguloACola(struct ColaRefinamiento *cola, struct Triangulo *t) {
    if(cola->numTriangulos < cola->capacidad) {
        cola->triangulos[cola->numTriangulos++] = t;
    }
}

struct Triangulo* extraerTriangulo(struct ColaRefinamiento *cola) {
    if(cola->numTriangulos > 0) {
        return cola->triangulos[--cola->numTriangulos];
    }
    return NULL;
}

bool dentroLimites(struct Punto *p, struct Triangulacion *tr) {
    // Encontrar los límites del dominio original
    double xmin = tr->puntos[0].x;
    double xmax = tr->puntos[0].x;
    double ymin = tr->puntos[0].y;
    double ymax = tr->puntos[0].y;
    
    for(int i = 1; i < tr->numPuntos; i++) {
        if(tr->puntos[i].x < xmin) xmin = tr->puntos[i].x;
        if(tr->puntos[i].x > xmax) xmax = tr->puntos[i].x;
        if(tr->puntos[i].y < ymin) ymin = tr->puntos[i].y;
        if(tr->puntos[i].y > ymax) ymax = tr->puntos[i].y;
    }
    
    // Margen más estricto
    const double margen = 0.05;
    
    // Verificación más estricta de los límites
    bool dentroX = p->x > xmin + margen && p->x < xmax - margen;
    bool dentroY = p->y > ymin + margen && p->y < ymax - margen;
    
    // Verificación adicional para puntos muy alejados
    bool noMuyLejos = fabs(p->y - ymin) < 0.2 && fabs(p->y - ymax) < 0.2;
    
    return dentroX && dentroY && noMuyLejos;
}

struct ListaTriangulos* encontrarTriangulosIntersectados(struct Triangulacion *tr, 
                                                        struct Punto *p1, 
                                                        struct Punto *p2) {
    // Inicializar lista de triángulos
    struct ListaTriangulos *lista = malloc(sizeof(struct ListaTriangulos));
    if (!lista) return NULL;
    
    lista->capacidad = 10;  // Capacidad inicial
    lista->numTriangulos = 0;
    lista->triangulos = malloc(lista->capacidad * sizeof(struct Triangulo*));
    if (!lista->triangulos) {
        free(lista);
        return NULL;
    }

    // Encontrar un triángulo inicial que contenga p1
    struct Triangulo *triangulo_inicial = NULL;
    for (int i = 0; i < tr->numTriangulos; i++) {
        if (puntoEnTriangulo(&tr->triangulos[i], p1)) {
            triangulo_inicial = &tr->triangulos[i];
            break;
        }
    }

    if (!triangulo_inicial) {
        free(lista->triangulos);
        free(lista);
        return NULL;
    }

    // Agregar el triángulo inicial a la lista
    lista->triangulos[lista->numTriangulos++] = triangulo_inicial;

    // Buscar triángulos intersectados
    int i = 0;
    while (i < lista->numTriangulos) {
        struct Triangulo *t = lista->triangulos[i];
        
        // Verificar vecinos
        for (int j = 0; j < 3; j++) {
            if (t->vecinos[j] != NULL) {
                // Verificar si el segmento intersecta con este triángulo
                double x, y;
                struct Punto *a = t->vertices[j];
                struct Punto *b = t->vertices[(j+1)%3];
                
                if (encontrarInterseccion(p1->x, p1->y, p2->x, p2->y,
                                        a->x, a->y, b->x, b->y,
                                        &x, &y)) {
                    // Verificar si ya está en la lista
                    bool ya_existe = false;
                    for (int k = 0; k < lista->numTriangulos; k++) {
                        if (lista->triangulos[k] == t->vecinos[j]) {
                            ya_existe = true;
                            break;
                        }
                    }
                    
                    // Agregar a la lista si no existe
                    if (!ya_existe) {
                        // Expandir la lista si es necesario
                        if (lista->numTriangulos >= lista->capacidad) {
                            lista->capacidad *= 2;
                            struct Triangulo **temp = realloc(lista->triangulos, 
                                                            lista->capacidad * sizeof(struct Triangulo*));
                            if (!temp) {
                                // Manejar error de memoria
                                free(lista->triangulos);
                                free(lista);
                                return NULL;
                            }
                            lista->triangulos = temp;
                        }
                        
                        lista->triangulos[lista->numTriangulos++] = t->vecinos[j];
                    }
                }
            }
        }
        i++;
    }

    printf("Encontrados %d triángulos intersectados\n", lista->numTriangulos);
    return lista;
}

bool compartenArista(struct Triangulo *t1, struct Triangulo *t2) {
    // Verificar si dos triángulos comparten una arista
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if ((t1->vertices[i] == t2->vertices[j] && 
                 t1->vertices[(i+1)%3] == t2->vertices[(j+1)%3]) ||
                (t1->vertices[i] == t2->vertices[(j+1)%3] && 
                 t1->vertices[(i+1)%3] == t2->vertices[j])) {
                return true;
            }
        }
    }
    return false;
}

void marcarAristasRestringidas(struct Triangulacion *tr, struct Punto *p1, struct Punto *p2) {
    // Marcar las aristas que forman parte del segmento como restringidas
    for (int i = 0; i < tr->numTriangulos; i++) {
        struct Triangulo *t = &tr->triangulos[i];
        for (int j = 0; j < 3; j++) {
            if ((t->vertices[j] == p1 && t->vertices[(j+1)%3] == p2) ||
                (t->vertices[j] == p2 && t->vertices[(j+1)%3] == p1)) {
                t->aristasRestringidas[j] = 1;
            }
        }
    }
}

struct Borde* actualizarBase(struct Triangulacion *tr, struct Borde *baseActual, struct Punto *nuevo) {
    struct Borde *nuevaBase = malloc(sizeof(struct Borde));
    if (!nuevaBase) return NULL;
    
    // La nueva base depende de qué lado se agregó el punto
    if (orientacion(baseActual->p1, baseActual->p2, nuevo) > 0) {
        // Punto agregado a la izquierda
        nuevaBase->p1 = nuevo;
        nuevaBase->p2 = baseActual->p2;
    } else {
        // Punto agregado a la derecha
        nuevaBase->p1 = baseActual->p1;
        nuevaBase->p2 = nuevo;
    }
    
    free(baseActual);
    return nuevaBase;
}

void crearSegmento(struct Triangulacion *tr, struct Punto *p1, struct Punto *p2) {
    if (tr->numBordes >= tr->maxBordes) return;
    
    struct Borde *borde = malloc(sizeof(struct Borde));
    if (!borde) return;
    
    borde->p1 = p1;
    borde->p2 = p2;
    tr->bordes[tr->numBordes++] = borde;
}


/* Funciones de triangulación                                                  */

// Función para inicializar la triangulación
struct Triangulacion* inicializarTriangulacion(struct Punto* puntos, int numPuntos, int numPuntosRegion1) {
    struct Triangulacion* tr = (struct Triangulacion*)malloc(sizeof(struct Triangulacion));
    if (!tr) return NULL;

    tr->maxTriangulos = numPuntos * 3;
    tr->triangulos = (struct Triangulo*)malloc(tr->maxTriangulos * sizeof(struct Triangulo));
    if (!tr->triangulos) {
        free(tr);
        return NULL;
    }

    tr->numTriangulos = 0;
    tr->puntos = puntos;
    tr->numPuntos = numPuntos;
    tr->maxPuntos = numPuntos;
    tr->numPuntosRegion1 = numPuntosRegion1;  // Guardamos el número de puntos de la región 1

    return tr;
}

// Función para crear un nuevo triángulo
struct Triangulo* crearTriangulo(struct Punto *v1, struct Punto *v2, struct Punto *v3) {
    struct Triangulo *t = malloc(sizeof(struct Triangulo));
    if (t == NULL) {
        printf("Error: No se pudo asignar memoria para el triángulo\n");
        return NULL;
    }
    
    // Permitir puntos artificiales (índices negativos) durante la triangulación
    t->vertices[0] = v1;
    t->vertices[1] = v2;
    t->vertices[2] = v3;
    t->indices[0] = v1->indice;
    t->indices[1] = v2->indice;
    t->indices[2] = v3->indice;
    
    // Marcar como triángulo artificial si tiene algún punto artificial
    if (v1->indice < 0 || v2->indice < 0 || v3->indice < 0) {
        t->esTrianguloSuper = 1;  // Usar 1 para indicar triángulo artificial
    } else {
        t->esTrianguloSuper = 0;
    }
    
    // Inicializar vecinos
    t->vecinos[0] = NULL;
    t->vecinos[1] = NULL;
    t->vecinos[2] = NULL;
    
    return t;
}

// Función para crear el super-triángulo que contendrá todos los puntos
struct Triangulo* crearSuperTriangulo(struct Punto *puntos, int numPuntos) {
    int i;
    // Encontrar los límites del conjunto de puntos
    double minX = puntos[0].x;
    double minY = puntos[0].y;
    double maxX = minX;
    double maxY = minY;
    
    for(i = 1; i < numPuntos; i++) {
        if(puntos[i].x < minX) minX = puntos[i].x;
        if(puntos[i].y < minY) minY = puntos[i].y;
        if(puntos[i].x > maxX) maxX = puntos[i].x;
        if(puntos[i].y > maxY) maxY = puntos[i].y;
    }
    
    // Calcular el centro y el tamaño del rectángulo que contiene los puntos
    double dx = maxX - minX;
    double dy = maxY - minY;
    double deltaMax = (dx > dy) ? dx : dy;
    double midX = (minX + maxX) / 2;
    double midY = (minY + maxY) / 2;

    // Crear los vértices del super-triángulo
    struct Punto *p1 = malloc(sizeof(struct Punto));
    struct Punto *p2 = malloc(sizeof(struct Punto));
    struct Punto *p3 = malloc(sizeof(struct Punto));
    
    // Hacer el super-triángulo lo suficientemente grande
    p1->x = midX - 20 * deltaMax;
    p1->y = midY - deltaMax;
    p1->indice = -1;
    
    p2->x = midX + 20 * deltaMax;
    p2->y = midY - deltaMax;
    p2->indice = -2;
    
    p3->x = midX;
    p3->y = midY + 20 * deltaMax;
    p3->indice = -3;
    
    // Crear el super-triángulo
    struct Triangulo *superTriangulo = crearTriangulo(p1, p2, p3);
    superTriangulo->esTrianguloSuper = 1;
    
    return superTriangulo;
}

void agregarTriangulo(struct Triangulacion* tr, struct Triangulo* t) {
    if (tr->numTriangulos >= tr->maxTriangulos) {
        // Aumentar capacidad si es necesario
        tr->maxTriangulos *= 2;
        struct Triangulo* temp = (struct Triangulo*)realloc(tr->triangulos, 
                                tr->maxTriangulos * sizeof(struct Triangulo));
        if (!temp) return;
        tr->triangulos = temp;
    }

    // Copiar el triángulo a la estructura
    tr->triangulos[tr->numTriangulos] = *t;
    tr->numTriangulos++;
}

void eliminarTriangulo(struct Triangulacion* tr, struct Triangulo* t) {
    // Buscar el triángulo en el arreglo
    int idx = -1;
    for (int i = 0; i < tr->numTriangulos; i++) {
        if (&tr->triangulos[i] == t) {
            idx = i;
            break;
        }
    }

    if (idx == -1) return;  // Triángulo no encontrado

    // Mover el último triángulo a la posición del eliminado
    if (idx < tr->numTriangulos - 1) {
        tr->triangulos[idx] = tr->triangulos[tr->numTriangulos - 1];
    }
    tr->numTriangulos--;
}

// Función para optimización local después de insertar un punto
void optimizarLocal(struct Triangulacion *tr, struct Punto *p) {
    bool cambios;
    do {
        cambios = false;
        // Encontrar todos los pares de triángulos adyacentes que comparten el punto p
        for(int i = 0; i < tr->numTriangulos; i++) {
            if(tieneVertice(&tr->triangulos[i], p)) {
                for(int j = 0; j < 3; j++) {
                    struct Triangulo *vecino = tr->triangulos[i].vecinos[j];
                    if(vecino != NULL && !vecino->esTrianguloSuper) {
                        if(necesitaIntercambio(&tr->triangulos[i], vecino)) {
                            intercambiarDiagonal(tr, &tr->triangulos[i], vecino);
                            cambios = true;
                        }
                    }
                }
            }
        }
    } while(cambios);
}

void dividirTriangulo(struct Triangulacion *tr, struct Triangulo *t, struct Punto *p) {
    struct Punto *a = t->vertices[0];
    struct Punto *b = t->vertices[1];
    struct Punto *c = t->vertices[2];

    // Crear tres nuevos triángulos
    struct Triangulo* t1 = crearTriangulo(p, a, b);
    struct Triangulo* t2 = crearTriangulo(p, b, c);
    struct Triangulo* t3 = crearTriangulo(p, c, a);

    // Preservar la información del triángulo super
    t1->esTrianguloSuper = t->esTrianguloSuper;
    t2->esTrianguloSuper = t->esTrianguloSuper;
    t3->esTrianguloSuper = t->esTrianguloSuper;

    // Agregar los nuevos triángulos
    agregarTriangulo(tr, t1);
    agregarTriangulo(tr, t2);
    agregarTriangulo(tr, t3);

    // Eliminar el triángulo original
    eliminarTriangulo(tr, t);
}

// Función para eliminar triángulos que contienen vértices del super-triángulo
void eliminarTriangulosSuper(struct Triangulacion *tr) {
    int j = 0;
    for (int i = 0; i < tr->numTriangulos; i++) {
        if (!tieneVerticeArtificial(&tr->triangulos[i])) {
            if (i != j) tr->triangulos[j] = tr->triangulos[i];
            j++;
        }
    }
    tr->numTriangulos = j;
}

// Función para liberar la memoria de la triangulación
void liberarTriangulacion(struct Triangulacion *tr) {
    if (tr != NULL) {
        if (tr->puntos != NULL) {
            free(tr->puntos);
        }
        if (tr->triangulos != NULL) {
            free(tr->triangulos);
        }
        if (tr->bordes != NULL) {
            for (int i = 0; i < tr->numBordes; i++) {
                free(tr->bordes[i]);
            }
            free(tr->bordes);
        }
        free(tr);
    }
}

// Función para insertar un segmento como restricción
void insertarSegmentoRestriccion(struct Triangulacion *tr, struct Segmento *seg) {
    struct Punto *p1 = &tr->puntos[seg->v1];
    struct Punto *p2 = &tr->puntos[seg->v2];
    
    // Encontrar triángulos que intersecta el segmento
    struct ListaTriangulos *triangulos = encontrarTriangulosIntersectados(tr, p1, p2);
    if (!triangulos) return;

    // Crear nueva arista restringida
    for (int i = 0; i < triangulos->numTriangulos - 1; i++) {
        struct Triangulo *t1 = triangulos->triangulos[i];
        struct Triangulo *t2 = triangulos->triangulos[i + 1];
        
        // Dividir triángulos si es necesario
        if (!compartenArista(t1, t2)) {
            struct Punto *p = calcularPuntoInterseccion(t1, t2, p1, p2);
            if (p) {
                dividirTriangulo(tr, t1, p);
                dividirTriangulo(tr, t2, p);
                optimizarLocal(tr, p);
            }
        }
    }

    // Marcar aristas como restringidas
    marcarAristasRestringidas(tr, p1, p2);
}

// Modificar la triangulación principal para incluir restricciones
struct Triangulacion* triangulacionDelaunayRestringida(struct EntradaPoly *entrada) {
    printf("Iniciando triangulación de Delaunay con restricciones...\n");
    
    // Paso 1: Triangulación inicial
    struct Triangulacion* tr = triangulacionDelaunay(entrada->vertices, 
                                                   entrada->numVertices,
                                                   entrada->numVertices);
    if (!tr) {
        printf("[ERROR] Falló la triangulación inicial\n");
        return NULL;
    }
    
    printf("Triangulación inicial completada.\n");
    printf("Insertando restricciones de segmentos...\n");
    
    // Paso 2: Insertar segmentos como restricciones
    for (int i = 0; i < entrada->numSegmentos; i++) {
        printf("Procesando segmento %d de %d\n", i + 1, entrada->numSegmentos);
        insertarSegmentoRestriccion(tr, &entrada->segmentos[i]);
    }
    
    printf("Restricciones insertadas.\n");
    printf("Actualizando relaciones de vecindad...\n");
    
    // Paso 3: Actualizar estructura final
    actualizarVecinos(tr);
    
    printf("Triangulación restringida completada.\n");
    return tr;
}

struct Borde* encontrarBaseInicial(struct Triangulacion *tr, int medio) {
    struct Borde *base = malloc(sizeof(struct Borde));
    if (!base) return NULL;
    
    // Encontrar el punto más a la derecha de la mitad izquierda
    struct Punto *p1 = &tr->puntos[medio];
    // Y el punto más a la izquierda de la mitad derecha
    struct Punto *p2 = &tr->puntos[medio + 1];
    
    base->p1 = p1;
    base->p2 = p2;
    return base;
}

struct Punto* encontrarCandidatoIzquierda(struct Triangulacion *tr, struct Borde *base) {
    struct Punto *mejor = NULL;
    double mejorAngulo = -DBL_MAX;
    
    // Buscar entre todos los puntos a la izquierda de la base
    for (int i = 0; i < tr->numPuntos; i++) {
        struct Punto *p = &tr->puntos[i];
        if (p == base->p1 || p == base->p2) continue;
        
        // Verificar si está a la izquierda de la base
        if (orientacion(base->p1, base->p2, p) > 0) {
            double angulo = calcularAngulo(base->p1, base->p2, p);
            if (angulo > mejorAngulo) {
                mejorAngulo = angulo;
                mejor = p;
            }
        }
    }
    
    return mejor;
}

struct Punto* encontrarCandidatoDerecha(struct Triangulacion *tr, struct Borde *base) {
    struct Punto *mejor = NULL;
    double mejorAngulo = -DBL_MAX;
    
    // Buscar entre todos los puntos a la derecha de la base
    for (int i = 0; i < tr->numPuntos; i++) {
        struct Punto *p = &tr->puntos[i];
        if (p == base->p1 || p == base->p2) continue;
        
        // Verificar si está a la derecha de la base
        if (orientacion(base->p1, base->p2, p) < 0) {
            double angulo = calcularAngulo(base->p1, base->p2, p);
            if (angulo > mejorAngulo) {
                mejorAngulo = angulo;
                mejor = p;
            }
        }
    }
    
    return mejor;
}

void triangular(struct Triangulacion *tr) {
    printf("Iniciando ordenamiento de puntos...\n");
    qsort(tr->puntos, tr->numPuntos, sizeof(struct Punto), compararPuntosX);
    printf("Puntos ordenados. Total puntos: %d\n", tr->numPuntos);
    
    printf("Iniciando división recursiva...\n");
    divideVencerasDelaunay(tr, 0, tr->numPuntos - 1);
    
    printf("Triangulación completada. Número de triángulos: %d\n", tr->numTriangulos);
}

void agregarTrianguloATriangulacion(struct Triangulacion *tr, struct Punto *v1, struct Punto *v2, struct Punto *v3) {
    // Verificar que no sea un triángulo degenerado
    double area = calcularAreaTriangulo2(v1, v2, v3);
    if (fabs(area) < 1e-10) {
        printf("Advertencia: Intento de agregar triángulo degenerado\n");
        return;
    }
    
    // Verificar que no exista ya el triángulo
    for (int i = 0; i < tr->numTriangulos; i++) {
        struct Triangulo *t = &tr->triangulos[i];
        if ((t->vertices[0] == v1 && t->vertices[1] == v2 && t->vertices[2] == v3) ||
            (t->vertices[0] == v2 && t->vertices[1] == v3 && t->vertices[2] == v1) ||
            (t->vertices[0] == v3 && t->vertices[1] == v1 && t->vertices[2] == v2)) {
            printf("Advertencia: Triángulo duplicado\n");
            return;
        }
    }
    
    // Agregar el triángulo
    if (tr->numTriangulos < tr->maxTriangulos) {
        tr->triangulos[tr->numTriangulos].vertices[0] = v1;
        tr->triangulos[tr->numTriangulos].vertices[1] = v2;
        tr->triangulos[tr->numTriangulos].vertices[2] = v3;
        tr->triangulos[tr->numTriangulos].esTrianguloSuper = false;
        tr->numTriangulos++;
        printf("Triángulo agregado: (%d, %d, %d)\n", 
               v1->indice, v2->indice, v3->indice);
    }
}

void combinarTriangulaciones(struct Triangulacion *tr, int inicio, int medio, int fin) {
    printf("Combinando submallas [%d-%d] y [%d-%d]\n", inicio, medio, medio+1, fin);
    
    struct Borde *baseInicial = encontrarBaseInicial(tr, medio);
    if (!baseInicial) {
        printf("ERROR: No se pudo encontrar la base inicial\n");
        return;
    }
    
    while (true) {
        struct Punto *candidatoIzq = encontrarCandidatoIzquierda(tr, baseInicial);
        struct Punto *candidatoDer = encontrarCandidatoDerecha(tr, baseInicial);
        
        if (!candidatoIzq && !candidatoDer) break;
        
        if (candidatoIzq && (!candidatoDer || esDelaunay(candidatoIzq, candidatoDer, baseInicial))) {
            agregarTrianguloATriangulacion(tr, candidatoIzq, baseInicial->p1, baseInicial->p2);
            baseInicial->p2 = candidatoIzq;
        } else if (candidatoDer) {
            agregarTrianguloATriangulacion(tr, candidatoDer, baseInicial->p1, baseInicial->p2);
            baseInicial->p1 = candidatoDer;
        }
    }
    
    free(baseInicial);
}

int compararPuntosX(const void *a, const void *b) {
    struct Punto *p1 = (struct Punto *)a;
    struct Punto *p2 = (struct Punto *)b;
    
    if (p1->x < p2->x) return -1;
    if (p1->x > p2->x) return 1;
    return 0;
}

void divideVencerasDelaunay(struct Triangulacion *tr, int inicio, int fin) {
    printf("Procesando segmento [%d, %d]\n", inicio, fin);
    
    // Caso base: 2 o 3 puntos
    if (fin - inicio + 1 <= 3) {
        if (fin - inicio + 1 == 3) {
            struct Punto *p1 = &tr->puntos[inicio];
            struct Punto *p2 = &tr->puntos[inicio + 1];
            struct Punto *p3 = &tr->puntos[inicio + 2];
            
            double area = calcularAreaTriangulo2(p1, p2, p3);
            if (fabs(area) > 1e-10 && !existeTriangulo(tr, p1, p2, p3)) {
                if (orientacion(p1, p2, p3) > 0) {
                    agregarTrianguloATriangulacion(tr, p1, p2, p3);
                } else {
                    agregarTrianguloATriangulacion(tr, p1, p3, p2);
                }
            }
        }
        return;
    }
    
    int medio = inicio + (fin - inicio) / 2;
    
    divideVencerasDelaunay(tr, inicio, medio);
    divideVencerasDelaunay(tr, medio + 1, fin);
    
    combinarTriangulacionesOptimizado(tr, inicio, medio, fin);
}

void combinarTriangulacionesOptimizado(struct Triangulacion *tr, int inicio, int medio, int fin) {
    printf("Combinando submallas optimizado [%d-%d] y [%d-%d]\n", inicio, medio, medio+1, fin);
    
    // Crear un array para marcar puntos ya procesados
    bool *procesado = calloc(tr->numPuntos, sizeof(bool));
    if (!procesado) {
        printf("Error: No se pudo asignar memoria para el array de procesados\n");
        return;
    }
    
    // Encontrar puntos extremos para la base inicial
    struct Punto *p_izq = NULL, *p_der = NULL;
    double x_min = INFINITY, x_max = -INFINITY;
    
    for (int i = inicio; i <= fin; i++) {
        if (tr->puntos[i].x < x_min) {
            x_min = tr->puntos[i].x;
            p_izq = &tr->puntos[i];
        }
        if (tr->puntos[i].x > x_max) {
            x_max = tr->puntos[i].x;
            p_der = &tr->puntos[i];
        }
    }
    
    if (!p_izq || !p_der) {
        free(procesado);
        return;
    }
    
    // Procesar puntos
    for (int i = inicio; i <= fin; i++) {
        if (!procesado[i]) {
            struct Punto *p_actual = &tr->puntos[i];
            
            // Buscar los dos puntos más cercanos ya procesados
            struct Punto *p1 = NULL, *p2 = NULL;
            double min_dist1 = INFINITY, min_dist2 = INFINITY;
            
            for (int j = inicio; j <= fin; j++) {
                if (i != j && procesado[j]) {
                    double dist = distanciaEntrePuntos(p_actual, &tr->puntos[j]);
                    if (dist < min_dist1) {
                        min_dist2 = min_dist1;
                        p2 = p1;
                        min_dist1 = dist;
                        p1 = &tr->puntos[j];
                    } else if (dist < min_dist2) {
                        min_dist2 = dist;
                        p2 = &tr->puntos[j];
                    }
                }
            }
            
            if (p1 && p2) {
                // Verificar orientación y criterio de Delaunay
                if (orientacion(p1, p2, p_actual) > 0 && 
                    !existeTriangulo(tr, p1, p2, p_actual)) {
                    agregarTrianguloATriangulacion(tr, p1, p2, p_actual);
                }
            }
            
            procesado[i] = true;
        }
    }
    
    free(procesado);
}

bool existeTriangulo(struct Triangulacion *tr, struct Punto *p1, struct Punto *p2, struct Punto *p3) {
    for (int i =0; i < tr->numTriangulos; i++) {
        struct Triangulo *t = &tr->triangulos[i];
        if ((t->vertices[0] == p1 && t->vertices[1] == p2 && t->vertices[2] == p3) ||
            (t->vertices[0] == p2 && t->vertices[1] == p3 && t->vertices[2] == p1) ||
            (t->vertices[0] == p3 && t->vertices[1] == p1 && t->vertices[2] == p2)) {
            return true;
        }
    }
    return false;
}

// Función principal que inicia el proceso
struct Triangulacion* triangulacionDelaunay(struct Punto *puntos, int numPuntos, int numPuntosRegion1) {
    printf("Iniciando triangulación de Delaunay...\n");
    
    // Ordenar puntos por coordenada x
    qsort(puntos, numPuntos, sizeof(struct Punto), compararPuntosX);
    
    // Inicializar la triangulación
    struct Triangulacion* tr = inicializarTriangulacion(puntos, numPuntos, numPuntosRegion1);
    if (!tr) {
        printf("Error al inicializar la triangulación\n");
        return NULL;
    }
    
    // Aplicar el algoritmo divide y vencerás
    divideVencerasDelaunay(tr, 0, numPuntos - 1);
    
    // Refinar la malla
    double anguloMinimo = 20.0 * M_PI / 180.0;  // 20 grados en radianes
    double areaMaxima = calcularAreaMaximaPermitida(tr);
    refinarMalla(tr, anguloMinimo, areaMaxima);
    
    printf("Triangulación completada.\n");
    return tr;
}


/* Funciones de refinamiento                                                   */

struct ColaRefinamiento* iniciarColaRefinamiento(int capacidad) {
    // Asignar memoria para la estructura principal
    struct ColaRefinamiento *cola = malloc(sizeof(struct ColaRefinamiento));
    if (cola == NULL) {
        printf("Error: No se pudo asignar memoria para la cola\n");
        return NULL;
    }
    
    // Asignar memoria para los arreglos de triángulos y bordes
    cola->triangulos = malloc(capacidad * sizeof(struct Triangulo*));
    cola->bordes = malloc(capacidad * sizeof(struct Borde*));
    
    // Verificar asignación de memoria
    if (cola->triangulos == NULL || cola->bordes == NULL) {
        printf("Error: No se pudo asignar memoria para los arreglos\n");
        free(cola);
        return NULL;
    }
    
    // Inicializar los contadores
    cola->numTriangulos = 0;  // Número actual de triángulos en la cola
    cola->numBordes = 0;      // Número actual de bordes en la cola
    cola->capacidad = capacidad;  // Capacidad máxima de la cola
    
    return cola;
}

void agregarPuntoATriangulacion(struct Triangulacion *tr, struct Punto *p) {
    // Si necesitamos más espacio
    if (tr->numPuntos >= tr->maxPuntos) {
        int nuevaCapacidad = tr->maxPuntos * 2;
        struct Punto *nuevosPuntos = realloc(tr->puntos, 
                                            nuevaCapacidad * sizeof(struct Punto));
        if (nuevosPuntos == NULL) {
            printf("Error: No se pudo expandir el arreglo de puntos\n");
            return;
        }
        tr->puntos = nuevosPuntos;
        tr->maxPuntos = nuevaCapacidad;
    }
    
    // Agregar el nuevo punto
    tr->puntos[tr->numPuntos] = *p;
    tr->puntos[tr->numPuntos].indice = tr->numPuntos;
    tr->numPuntos++;
    
    printf("Nuevo punto agregado en (%f, %f), índice %d\n", 
           p->x, p->y, tr->numPuntos - 1);
}

void liberarColaRefinamiento(struct ColaRefinamiento *cola) {
    if(cola != NULL) {
        free(cola->triangulos);
        free(cola->bordes);
        free(cola);
    }
}

bool estaDentroDeLimites(struct Triangulacion *tr, struct Punto *p) {
    int dentro = 0;
    for (int i = 0; i < tr->numBordes; i++) {
        struct Borde *b = tr->bordes[i];
        if ((b->p1->y > p->y) != (b->p2->y > p->y) &&
            p->x < (b->p2->x - b->p1->x) * (p->y - b->p1->y) / 
                   (b->p2->y - b->p1->y) + b->p1->x) {
            dentro = !dentro;
        }
    }
    return dentro;
}

bool hayPuntoCercano(struct Triangulacion *tr, struct Punto *p) {
    const double DISTANCIA_MINIMA = 0.001;
    
    for (int i = 0; i < tr->numPuntos; i++) {
        double dx = tr->puntos[i].x - p->x;
        double dy = tr->puntos[i].y - p->y;
        double dist = sqrt(dx*dx + dy*dy);
        
        if (dist < DISTANCIA_MINIMA) {
            return true;
        }
    }
    return false;
}

double calcularAreaMaximaPermitida(struct Triangulacion *tr) {
    double xmin = tr->puntos[0].x;
    double xmax = tr->puntos[0].x;
    double ymin = tr->puntos[0].y;
    double ymax = tr->puntos[0].y;
    
    for (int i = 1; i < tr->numPuntos; i++) {
        if (tr->puntos[i].x < xmin) xmin = tr->puntos[i].x;
        if (tr->puntos[i].x > xmax) xmax = tr->puntos[i].x;
        if (tr->puntos[i].y < ymin) ymin = tr->puntos[i].y;
        if (tr->puntos[i].y > ymax) ymax = tr->puntos[i].y;
    }
    
    double areaTotal = (xmax - xmin) * (ymax - ymin);
    return areaTotal * 0.01; // 1% del área total
}

void refinarMalla(struct Triangulacion *tr, double anguloMinimo, double areaMaxima) {
    printf("\n=== INICIO DEL REFINAMIENTO ===\n");
    int puntosIniciales = tr->numPuntos;
    int iteraciones = 0;
    const int MAX_ITERACIONES = 100;
    bool seAgregaronPuntos;
    
    do {
        seAgregaronPuntos = false;
        iteraciones++;
        
        // Lista temporal de puntos a agregar
        struct Punto *nuevosPuntos = malloc(tr->numTriangulos * sizeof(struct Punto));
        int numNuevosPuntos = 0;
        
        // Revisar cada triángulo
        for (int i = 0; i < tr->numTriangulos; i++) {
            struct Triangulo *t = &tr->triangulos[i];
            if (t->esTrianguloSuper) continue;
            
            // Calcular ángulos del triángulo
            double angulos[3];
            calcularAngulos(t, angulos);
            
            // Verificar si necesita refinamiento
            bool necesitaRefinar = false;
            for (int j = 0; j < 3; j++) {
                if (angulos[j] < anguloMinimo) {
                    necesitaRefinar = true;
                    break;
                }
            }
            
            if (necesitaRefinar || calcularAreaTriangulo(t) > areaMaxima) {
                struct Punto *circuncentro = calcularCircuncentro(t);
                if (circuncentro && estaDentroDeLimites(tr, circuncentro) && 
                    !hayPuntoCercano(tr, circuncentro)) {
                    nuevosPuntos[numNuevosPuntos++] = *circuncentro;
                    seAgregaronPuntos = true;
                }
                free(circuncentro);
            }
        }
        
        // Agregar los nuevos puntos y retriangular
        if (seAgregaronPuntos) {
            for (int i =0; i < numNuevosPuntos; i++) {
                if (tr->numPuntos < tr->maxPuntos) {
                    tr->puntos[tr->numPuntos++] = nuevosPuntos[i];
                    printf("Punto Steiner agregado: (%f, %f)\n", 
                           nuevosPuntos[i].x, nuevosPuntos[i].y);
                }
            }
            
            // Retriangular desde cero con los nuevos puntos
            divideVencerasDelaunay(tr, 0, tr->numPuntos - 1);
        }
        
        free(nuevosPuntos);
        
    } while (seAgregaronPuntos && iteraciones < MAX_ITERACIONES);
    
    printf("Refinamiento completado: %d puntos agregados en %d iteraciones\n",
           tr->numPuntos - puntosIniciales, iteraciones);
}

void imprimirEstadisticas(struct Triangulacion *tr) {
    printf("\nEstadísticas de la triangulación:\n");
    printf("Número de puntos: %d\n", tr->numPuntos);
    printf("Número de triángulos: %d\n", tr->numTriangulos);
    printf("Número de bordes: %d\n", tr->numBordes);
}

void verificarTriangulacion(struct Triangulacion *tr) {
    printf("\nVerificando triangulación...\n");
    int triangulos_invalidos = 0;
    
    for (int i = 0; i < tr->numTriangulos; i++) {
        // Verificar que los vértices sean distintos
        if ((tr->triangulos[i].vertices[0] == tr->triangulos[i].vertices[1] ||
            tr->triangulos[i].vertices[1] == tr->triangulos[i].vertices[2] ||
            tr->triangulos[i].vertices[2] == tr->triangulos[i].vertices[0])) {
            triangulos_invalidos++;
        }
    }
    
    printf("Triángulos inválidos encontrados: %d\n", triangulos_invalidos);
}


/* Funciones de interacción con el usuario                                     */
void info(){
    printf("\nInformacion sobre el programa:\n");
    printf("    -p  Triangula un grafo planar de lineas (.poly file).\n");
    printf("    -r  Refina una malla previamente generada.\n");
    printf("    -q  Genera una malla de calidad. Se puede especificar un angulo minimo.\n");
    printf("    -a  Aplica una restriccion de area maxima a los triangulos.\n");
    printf("    -D  Conforme a Delaunay: todos los triangulos son verdaderamente Delaunay.\n");
    printf("    -i: Muestra esta informacion\n");
    printf("    -s: Salir del programa\n");
    printf("\nPresione Enter para continuar...");
    getchar();
    getchar();
    system("cls");
    return;
}

void menu(){
    printf("==================================================\n");
    printf("             TRIANGULACION DE DELAUNAY            \n");
    printf("==================================================\n");
    printf("Este programa genera una triangulacion de Delaunay\n");
    printf("==================================================\n");
    printf("Seleccione una opcion:\n");
    printf("\n");
    printf("-p [archivo.poly]: Generar triangulacion\n");
    printf("-r: Refinar malla\n");
    printf("-i: Mostrar informacion\n");
    printf("-s: Salir\n");
    printf("==================================================\n");
    printf("Ingrese un comando: ");
}

// Agregar estas implementaciones
double calcularAreaTriangulo2(struct Punto *p1, struct Punto *p2, struct Punto *p3) {
    return fabs((p2->x - p1->x) * (p3->y - p1->y) - 
                (p3->x - p1->x) * (p2->y - p1->y)) / 2.0;
}

struct Borde* crearBorde(struct Triangulacion *tr, struct Punto *p1, struct Punto *p2) {
    struct Borde *borde = malloc(sizeof(struct Borde));
    if (!borde) return NULL;
    
    borde->p1 = p1;
    borde->p2 = p2;
    return borde;
}

void agregarBorde(struct Triangulacion *tr, struct Punto *p1, struct Punto *p2) {
    if (tr->numBordes >= tr->maxBordes) {
        int nuevaCapacidad = tr->maxBordes * 2;
        struct Borde **nuevosBordes = realloc(tr->bordes, 
                                             nuevaCapacidad * sizeof(struct Borde*));
        if (!nuevosBordes) return;
        
        tr->bordes = nuevosBordes;
        tr->maxBordes = nuevaCapacidad;
    }
    
    struct Borde *nuevoBorde = crearBorde(tr, p1, p2);
    if (nuevoBorde) {
        tr->bordes[tr->numBordes++] = nuevoBorde;
    }
}

int main() {
    char comando[100];
    char nombreArchivo[100];
    struct Triangulacion *tr = NULL;
    struct EntradaPoly *entrada = NULL;
    
    do {
        menu();
        scanf("%s", comando);
        
        if (strcmp(comando, "-p") == 0) {
            printf("Ingrese el nombre del archivo .poly: ");
            scanf("%s", nombreArchivo);
            
            // Leer archivo .poly
            entrada = leerArchivoPoly(nombreArchivo);
            if (!entrada) {
                printf("[ERROR] No se pudo procesar el archivo %s\n", nombreArchivo);
                continue;
            }

            // Inicializar la triangulación
            tr = malloc(sizeof(struct Triangulacion));
            if (!tr) {
                printf("[ERROR] No se pudo asignar memoria para la triangulación\n");
                liberarEntradaPoly(entrada);
                continue;
            }

            // Configurar la triangulación
            tr->maxPuntos = entrada->numVertices;
            tr->numPuntos = entrada->numVertices;
            tr->puntos = malloc(tr->maxPuntos * sizeof(struct Punto));
            tr->maxTriangulos = 2 * tr->maxPuntos;
            tr->triangulos = malloc(tr->maxTriangulos * sizeof(struct Triangulo));
            tr->numTriangulos = 0;
            tr->maxBordes = 3 * tr->maxPuntos;
            tr->bordes = malloc(tr->maxBordes * sizeof(struct Borde*));
            tr->numBordes = 0;

            if (!tr->puntos || !tr->triangulos || !tr->bordes) {
                printf("[ERROR] No se pudo asignar memoria para las estructuras\n");
                liberarTriangulacion(tr);
                liberarEntradaPoly(entrada);
                continue;
            }

            // Copiar puntos del archivo .poly
            for (int i = 0; i < entrada->numVertices; i++) {
                tr->puntos[i].x = entrada->vertices[i].x;
                tr->puntos[i].y = entrada->vertices[i].y;
                tr->puntos[i].indice = i;
            }

            // Realizar la triangulación usando divide y vencerás
            printf("\nIniciando triangulación de Delaunay...\n");
            triangular(tr);  // Esta función ya ordena los puntos y aplica divide y vencerás
            printf("Triangulación básica completada.\n");

            // Agregar restricciones de bordes
            printf("\nAgregando restricciones de bordes...\n");
            for (int i = 0; i < entrada->numSegmentos; i++) {
                insertarSegmentoRestriccion(tr, &entrada->segmentos[i]);
            }
            printf("Restricciones de bordes completadas.\n");

            // Actualizar estructura final
            actualizarVecinos(tr);
            
            imprimirEstadisticas(tr);
            verificarTriangulacion(tr);

            // Generar archivos de salida
            char nombreSalida[256];
            strcpy(nombreSalida, nombreArchivo);
            nombreSalida[strlen(nombreSalida) - 5] = '\0'; // Quitar '.poly'
            
            char archivoNode[256], archivoEle[256];
            sprintf(archivoNode, "%s.node", nombreSalida);
            sprintf(archivoEle, "%s.ele", nombreSalida);
            
            guardarArchivoNode(tr, archivoNode);
            guardarArchivoEle(tr, archivoEle);

            printf("\nArchivos generados exitosamente:\n");
            printf("- %s\n", archivoNode);
            printf("- %s\n", archivoEle);

            // Liberar memoria
            liberarTriangulacion(tr);
            liberarEntradaPoly(entrada);
            
            printf("\nPresione Enter para continuar...");
            getchar();
            getchar();
        }
        else if (strcmp(comando, "-r") == 0) {
            printf("\nRefinando malla...\n");
            // TODO: Implementar refinamiento
            printf("\nPresione Enter para continuar...");
            getchar();
            getchar();
        }
        else if (strcmp(comando, "-i") == 0) {
            info();
        }
        else if (strcmp(comando, "-s") == 0) {
            printf("\nSaliendo del programa...\n");
            if (tr) liberarTriangulacion(tr);
            if (entrada) liberarEntradaPoly(entrada);
            break;
        }
        else {
            printf("\nComando no válido. Por favor, intente nuevamente.\n");
            printf("\nPresione Enter para continuar...");
            getchar();
            getchar();
        }
        
        system("cls");  // Limpiar pantalla

    } while(1);
    
    return 0;
}