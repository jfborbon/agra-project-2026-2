/*
Proyecto Final AGRA 2026-2

Integrantes:
    - Juan Felipe Borbón Melo
    - Juan Andrés Correa Arenas

Análisis de complejidad:
... TODO 
*/
#include <vector>
#include <queue>
#include <climits>
#include <stdio.h>

using namespace std;

// abreviación para la representación de grafo ponderado por medio de lista de aristas
typedef vector<vector<pair<int, int>>> Graph;

/*
Dijkstra que se usa para calcular las distancias mínimas hacia cada uno de los nodos,
este se termina empleando tanto desde la posición de Lion-O como la de Mumm-Ra. Para
el caso de Lion-O debido a que independientemente de la velocidad que este pueda llegar
a alcanzar este no sera capaz de pasar por una intersección en la cual esta Mumm-Ra
entonces se usa el entero `blocked` y la condición `if (nextInter != blocked)` para
unicamente procesar las aristas que atraviesen dicha intersección bloqueada. Para el caso
de Mumm-Ra simplemente se bloquea una intersección inexistente con algún valor no
alcanzable (como -1).
*/
vector<int> dijkstra(Graph& graph, int start, int blocked) {
    
    vector<int> dist(graph.size(), INT_MAX);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    dist[start] = 0;
    pq.push({0, start});

    while (!pq.empty()) {
        int intersection = pq.top().second;
        int intersectionCost = pq.top().first;
        pq.pop();
        if (intersectionCost == dist[intersection]) {
            for (pair<int, int>& nextIntersection : graph[intersection]) {
                int nextInter = nextIntersection.first;
                if (nextInter != blocked) { // en caso de que esa intersección este bloqueada (como el de Lion-O para la posición de Mumm-Ra) esta intersección no se procesa
                    int nextInterDistance = nextIntersection.second;
                    int cost = nextInterDistance + intersectionCost;
                    if (cost < dist[nextInter]) {
                        pq.push({cost, nextInter});
                        dist[nextInter] = cost;
                    }
                }
            }
        }
    }

    return dist;
}


/*
Función para encontrar el limite máximo para la velocidad de Lion-O, se busca encontrar la
velocidad minima para la cual lion-O llega al mismo tiempo que Mumm-Ra para cualquiera de las
intersecciones y se toma la mayor de estas. Esto se basa en que la velocidad que Lion-O necesita
para llegar a la par que Mumm-Ra es

Velocidad Lion-O = distancia Lion-O / tiempo Mumm-Ra
*/
long double lionOMaxSpeed(vector<int>& lionODistances, vector<long double>& mummRaTimes, int mummRa) {
    long double ans = -1.0L;
    for (int intersection = 0; intersection < lionODistances.size(); intersection++) {
        if (intersection != mummRa && lionODistances[intersection] != INT_MAX) {
            ans = max(ans, (long double)lionODistances[intersection] / mummRaTimes[intersection]);
        }
    }
    return ans;
}


/*
Función que permite determinar si Lion-O es capaz de escapar de Mumm-Ra con una velocidad
dada para Lion-O. Es un Dijkstra que parte desde la posición inicial de Lion-O y que solo
encola un nodo si Lion-O puede llegar en menor tiempo a ese lugar. De esta forma si Lion-O
logra llegar a uno de los puntos de escape primero entonces se da por hecho que es posible
escapar y se retorna `true`. Caso contrario si Lion-O no pudo encontrar ningún camino que
le permitiera llegar a la meta sin ser atrapado entonces se retorna `false`.
*/
bool lionOCanScape(Graph& graph, vector<bool>& isScapePoint, vector<long double>& mummRaTimes, int lionO, long double lionOSpeed) {
    bool ans = false;
    vector<int> dist(graph.size(), INT_MAX);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    dist[lionO] = 0;
    pq.push({0, lionO});

    while (!pq.empty() && !ans) {
        int intersection = pq.top().second;
        int intersectionCost = pq.top().first;
        pq.pop();

        if (isScapePoint[intersection]) { // si la intersección es un punto de escape se da por echo que para esta velocidad Lion-O si puede escapar
            ans = true;
        }
        else if (intersectionCost == dist[intersection]) {
            for (pair<int, int>& nextIntersection : graph[intersection]) {
                int nextInter = nextIntersection.first;
                int nextInterDistance = nextIntersection.second;
                int cost = nextInterDistance + intersectionCost;
                if (cost < dist[nextInter] && (long double)cost / lionOSpeed <= mummRaTimes[nextInter]) { // se procesa la intersección si Lion-O puede llegar primero que Mumm-Ra
                    pq.push({cost, nextInter});
                    dist[nextInter] = cost;
                }
            }
        }
    }

    return ans;
}


/*
Función que se encarga de realizar bisección sobre la velocidad requerida para que Lion-O
pueda escapar, en otras palabras va tanteando con que velocidades puede escapar, si ve que
con una velocidad dada Lion-O puede escapar entonces prueba con una que sea menor, por el
contrario si para una velocidad dada encuentra que Lion-O no puede escapar prueba con una
que sea mayor. Todo esto gracias a los limites left y right para aplicar una bisección e ir
reduciendo el rango de búsqueda cada vez a la mitad. 
*/
long double lionOScapeSpeedAux(
    Graph& graph, 
    vector<bool>& isScapePoint, 
    int lionO,
    vector<long double>& mummRaTimes,
    long double left,
    long double right, 
    long double error) {
    
    long double ans;
    if (right - left <= error) {
        ans = right;
    }
    else {
        long double mid = (right + left) / 2.0L;
        if (lionOCanScape(graph, isScapePoint, mummRaTimes, lionO, mid)) {
            right = mid;
        }
        else {
            left = mid;
        }
        ans = lionOScapeSpeedAux(graph, isScapePoint, lionO, mummRaTimes, left, right, error);
    }

    return ans;
}


/*
Función que se encarga de determinar si existe un camino posible para que Lion-O pueda 
llegar a uno de los puntos de salida. Para luego sacar el limite máximo y aplicar la 
bisección sobre las posibles velocidades que podría usar Lion-O para escapar de Mumm-Ra.
*/
long double lionOScapeSpeed(
    Graph& graph,
    vector<bool>& isScapePoint,
    vector<int>& scapePoints,
    int lionO,
    int mummRa,
    int mummRaSpeed) {
    
    long double ans = -1.0L; // Se inicializa el valor de retorno en -1 en cazo de que no se encuentre una salida para Lion-O
    vector<int> lionODistances = dijkstra(graph, lionO, mummRa);

    // Comprueba si existe algún camino de salida previamente antes de
    // realizar el procedimiento evaluando si es posible llegar a al
    // menos un punto de escape.
    size_t i = 0;
    bool thersAScape = false;
    while (i < scapePoints.size() && !thersAScape) {
        if (lionODistances[scapePoints[i]] != INT_MAX) {
            thersAScape = true;
        }
        i++;
    }
    if (thersAScape) {
        vector<int> mummRaDistances = dijkstra(graph, mummRa, -1);

        long double error = 0.0000001L; // se usa un error de 10^{-7}

        // se saca el tiempo que toma mummRa para llegar a cualquier intersección
        vector<long double> mummRaTimes;
        mummRaTimes.reserve(graph.size());
        for (int mummRaDist : mummRaDistances) {
            mummRaTimes.push_back((long double)mummRaDist / (long double)mummRaSpeed);
        }

        long double maxSpeed = lionOMaxSpeed(lionODistances, mummRaTimes, mummRa);
        
        ans = lionOScapeSpeedAux(graph, isScapePoint, lionO, mummRaTimes, 0.0L, maxSpeed, error);
    }

    return ans;
}


int main() {

    int T, n, m, k, v, a, b, d, s, r, ki;
    long double result;

    scanf("%d", &T);
    while (T > 0) {

        scanf("%d %d %d %d", &n, &m, &k, &v);

        Graph graph(n);
        while (m > 0) {
            scanf("%d %d %d", &a, &b, &d);
            a--; b--;
            graph[a].push_back({b, d});
            graph[b].push_back({a, d});
            m--;
        }

        scanf("%d %d", &s, &r);
        s--; r--;
        vector<bool> isScapePoint(n, false);
        vector<int> scapePoints;
        scapePoints.reserve(k);
        while (k > 0) {
            scanf("%d", &ki);
            ki--;
            isScapePoint[ki] = true;
            scapePoints.push_back(ki);
            k--;
        }

        result = lionOScapeSpeed(graph, isScapePoint, scapePoints, s, r, v);

        if (result != -1.0L) {
            printf("Lion-O can escape at a speed of %.7Lf km/h!\n", result);
        }
        else {
            printf("Mumm-Ra captures Lion-O!\n");
        }

        T--;
    }

    return 0;
}