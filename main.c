#include <stdio.h>
#include "stdbool.h"
#include "stdlib.h"
#include "string.h"
#include <math.h>
#include <limits.h>
#include <time.h>
#define  MAXROW 10
#define  MAXLINE 10

#define INIT_SIZE 10
#define MULTIPLIER 2
#define STACK_OVERFLOW  -100
#define STACK_UNDERFLOW -101
#define OUT_OF_MEMORY   -102

// structure for paka
typedef struct Vector {
    int y;
    int *vector;
} Vector;


typedef struct Maze {
    int size;
    int **maze;
} Maze;


typedef struct _Point {
    int _x; // map point abscissa
    int _y; // Map point ordinate
} Point;

typedef struct _Node // queue node
{
    Point data;
    struct _Node *next;

} Node;

typedef struct _Queue // Queue
{
    Node *front;
    Node *rear;
} Queue;

typedef struct Result {
    int size;
    int count;
    Point *way;
} Result;


void InitQueue(Queue *q); // initialization
bool IsQueueEmpty(Queue *q); // appreciate the emptiness
void enQueue(Queue *q, Point ch); // place in queue
Point deQueue(Queue *q); // Dequeuing
void ClearQueue(Queue *q); // clear



void InitQueue(Queue *q) // initialization
{
    q->front = q->rear = (Node *) malloc(sizeof(Node));
    q->front->next = q->rear->next = NULL;
}

bool IsQueueEmpty(Queue *q)
{
    return q->front == q->rear;
}


void enQueue(Queue *q, Point ch) // place in queue
{
    Node *tmp = (Node *) malloc(sizeof(Node));
    tmp->data = ch;
    q->rear->next = tmp;
    q->rear = tmp;
    tmp->next = NULL;
}

Point deQueue(Queue *q) // Dequeuing
{
    Point data = q->front->next->data;
    if (q->front->next == q->rear) {
        q->rear = q->front;
        free(q->front->next);
        q->front->next = NULL;
    } else {
        Node *tmp = q->front->next;
        q->front->next = tmp->next;
        free(tmp);
    }
    return data;
}

void ClearQueue(Queue *q)
{
    Node *head = q->front->next;
    q->rear->next = NULL;
    q->rear = q->front;
    Node *t = NULL;
    while (head) {
        t = head->next;
        free(head);
        head = t;
    }
//    free(q->rear);
//    free(q->front);
}

void Show(Queue *q) {
    Node *tmp = q->front->next;
    while (tmp) {
        printf("(%d,%d)\t", tmp->data._x, tmp->data._y);
        tmp = tmp->next;
    }
    printf("\n");
}


Point prePoint [MAXROW] [MAXLINE]; // Used to store 2D effective path coordinates


void DisplayPre() // Printing the coordinates of the distance traveled
{
    for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 7; j++) {
            printf("(%2d,%2d)", prePoint[i][j]);
        }
        printf("\n");
    }
}

typedef int T;
typedef struct Stack_tag {
    T *DATA;
    size_t size;
    size_t top;
} Stack_t;

Stack_t *createStack() {
    Stack_t *out = NULL;
    out = (Stack_t *) malloc(sizeof(Stack_t));
    if (out == NULL) {
        exit(OUT_OF_MEMORY);
    }
    out->size = INIT_SIZE;
    out->DATA = (int *) malloc(out->size * sizeof(T));
    if (out->DATA == NULL) {
        free(out);
        exit(OUT_OF_MEMORY);
    }
    out->top = 0;
    return out;
}

void deleteStack(Stack_t **stack) {
    free((*stack)->DATA);
    free(*stack);
    *stack = NULL;
}

void resize(Stack_t *stack) {
    stack->size *= MULTIPLIER;
    stack->DATA = (int *) realloc(stack->DATA, stack->size * sizeof(T));
    if (stack->DATA == NULL) {
        exit(STACK_OVERFLOW);
    }
}

void push(Stack_t *stack, T value) {
    if (stack->top >= stack->size) {
        resize(stack);
    }
    stack->DATA[stack->top] = value;
    stack->top++;
}

T pop(Stack_t *stack) {
    if (stack->top == 0) {
        exit(STACK_UNDERFLOW);
    }
    stack->top--;
    return stack->DATA[stack->top];
}

T peek(const Stack_t *stack) {
    if (stack->top <= 0) {
        exit(STACK_UNDERFLOW);
    }
    return stack->DATA[stack->top - 1];
}

void implode(Stack_t *stack) {
    stack->size = stack->top;
    stack->DATA = (int *) realloc(stack->DATA, stack->size * sizeof(T));
}


void printMaze(int **maze, int n) // Print the map of the labyrinth
{
    for (int i = 0; i < n + 2; i++) {
        for (int j = 0; j < n + 2; j++) {
            if (maze[i][j] == 1) // 1 means wall
                printf("%2s", "*");
            else if (maze[i][j] == 2)
                printf("%2s", "#");
            else if (maze[i][j] == 3)
                printf("%2s", "="); // # Indicates the path to be taken
            else // 2 means the distance traveled
                printf("%2s", " ");
        }
        printf("\n");
    }
    printf("\n");

}

int g = 0;

int **copy_maze(int **maze, int n) {
    int **new_maze = (int **) malloc(sizeof(int *) * (n + 2));
    for (int i = 0; i < n + 2; i++) {
        new_maze[i] = (int *) malloc(sizeof(int) * (n + 2));
        for (int j = 0; j < n + 2; j++) {
            new_maze[i][j] = maze[i][j];
        }
    }
    return new_maze;
}

// Assign the coordinates of the points corresponding to the path to the queue, and assign the last pushed coordinates to the two-dimensional array
void visited(int i, int j, Point t, Queue *q, Point **prePoint) {
    g++;
    Point tmp;
    tmp._x = j;
    tmp._y = i;
    enQueue(q, tmp);
    prePoint[j][i] = t;
}


void clear_maze(int **maze, int n) {
    for (int i = 0; i < n + 2; ++i) {
        free(maze[i]);
    }
    free(maze);
}

int **init_maze(int n) {
    int **maze = (int **) malloc(sizeof(int *) * (n + 2));
    for (int i = 0; i < n + 2; i++) {
        maze[i] = (int *) malloc(sizeof(int) * (n + 2));
        for (int j = 0; j < n + 2; j++) {
            char d;
            if (i == 0 || j == 0 || i == n + 1 || j == n + 1) {
                maze[i][j] = 1;
            } else {
                scanf(" %c", &d);
                maze[i][j] = d - 48;
            }
        }
    }
    int **copy = copy_maze(maze, n);
    for (int i = 0; i < n + 2; ++i) {
        for (int j = 0; j < n + 2; ++j) {
            maze[i][j] = copy[j][i];
        }
    }
    clear_maze(copy, n);
    return maze;
}


int **init_maze_without_scanf(int n) {
    int **maze = (int **) malloc(sizeof(int *) * (n + 2));
    for (int i = 0; i < n + 2; i++) {
        maze[i] = (int *) malloc(sizeof(int) * (n + 2));
        for (int j = 0; j < n + 2; j++) {
            if (i == 0 || j == 0 || i == n + 1 || j == n + 1) {
                maze[i][j] = 1;
            } else {
                maze[i][j] = 0;
            }
        }
    }
    return maze;
}


void clear_ways(int **maze, int n) {
    for (int l = 0; l < n + 2; ++l) {
        for (int m = 0; m < n + 2; ++m) {
            maze[l][m] = maze[l][m] == 2 ? 0 : maze[l][m];
        }
    }
}


Result check_maze(int **maze, int n, int position_pak, Point EndPoint) {
    // Used to store 2D effective path coordinates
    Point **prePoint = (Point **) malloc(sizeof(Point *) * (n + 2));
    for (int i = 0; i < n + 2; ++i) {
        prePoint[i] = (Point *) malloc(sizeof(Point) * (n + 2));
        for (int j = 0; j < n + 2; ++j) {
            prePoint[i][j] = (Point) {-1, -1};
        }
    }
    Queue q; // Defining the global stack queue
    InitQueue(&q); // Initialize the queue
    Point StartPoint = {1, 1}; // Координаты входа в лабиринт
    if (maze[1][1] == 1) {
        Result r;
        r.size = -1;
        r.count = 0;
        return r;
    }
    enQueue(&q, StartPoint); // We put the initial coordinates in the queue
    int flag = 0; // Set a flag bit to determine if the exit is finally found
    int k = 0;
    while (!IsQueueEmpty(&q)) {
        k++;
        Point t = deQueue(&q);
        maze[t._x][t._y] = 2; // The queue coordinates are the points that can be traversed

//        printf("%d=%d %d=%d\n", t._x, EndPoint._x, t._y, EndPoint._y);
        // If queue coordinates are maze exit coordinates
        if (t._x == EndPoint._x && t._y == EndPoint._y) {
            flag = 1; // Table flag is 0 to mark exit
            break;
        }
//        printMaze(maze, n);
        // Determine if the right side of the queue coordinates is a valid node
        if (t._x + 1 <= n + 2 && maze[t._x + 1][t._y] == 0) {
            visited(t._x + 1, t._y, t, &q, prePoint);
            maze[t._x + 1][t._y] = 3;
        }
        // Determine if left side of queue coordinates is a valid node
        if (t._x - 1 >= 0 && maze[t._x - 1][t._y] == 0) {
            visited(t._x - 1, t._y, t, &q, prePoint);
            maze[t._x - 1][t._y] = 3;
        }
        // Determine if top of queue coordinates is a valid node
        if (t._y + 1 <= n + 2 && maze[t._x][t._y + 1] == 0) {
            visited(t._x, t._y + 1, t, &q, prePoint);
            maze[t._x][t._y + 1] = 3;
        }
        // Determine if bottom of queue coordinates is a valid node
        if (t._y - 1 >= 0 && maze[t._x][t._y - 1] == 0) {
            visited(t._x, t._y - 1, t, &q, prePoint);
            maze[t._x][t._y - 1] = 3;
        }

    }
    ClearQueue(&q); // Clearing the stack area

    if (flag == 0) {
        Result r;
        r.size = -1;
        r.count = 0;
        clear_ways(maze, n);
        return r;
    }
//    DisplayPre();
    Stack_t *s = createStack();
    Point t = EndPoint;
    int count = -1;
    while (t._x > 0) // Point forward coordinate until input -1 is found
    {
        count++;
        push(s, t._y);
        push(s, t._x);
        t = prePoint[t._x][t._y];
    }

    implode(s);
    Point *way = (Point *) malloc(sizeof(Point) * count);
    for (int i = 0; i < count + 1; i++) {
        int Y = pop(s);
        int X = pop(s);
        Point p = {X, Y};
        way[i] = p;
    }

    deleteStack(&s);
    Result r;
    r.size = count + position_pak * 2;
    r.way = way;
    r.count = count;
    clear_ways(maze, n);
    return r;
}


int **apply_pak(Vector paka, int **maze, int n) {
    int **maze_after_pak = init_maze_without_scanf(n);
    for (int i = 1; i < n + 1; i++) {
        for (int j = 1; j < n + 1; j++) {
            if (paka.vector[j - 1] == 1) {
                maze_after_pak[i][j] = maze[i][j] == 1 ? 0 : 1;
            } else {
                maze_after_pak[i][j] = maze[i][j];
            }
        }
    }
    clear_maze(maze, n);
    return maze_after_pak;
}

int **next_combination_pak(Vector paka, int n) {
//    10pak! = 3628800 combination
//    new_combination =
//        return new_combination;
}


int *convert_base(long value, int base, int *output_number) { // 21
    int *array = (int *) calloc(sizeof(int), base); // [0,?,?]

    while (value >= base) {
        array[*output_number] = (int) (value % base);
        value = value / base;
        (*output_number)++;
    }

    array[*output_number] = (int) value;
    (*output_number)++;
    return array;
}

// [0417520] => [0 => 2, 1 => 1, 2 => 1, 3 => 0, 4 => 1, 5=>1, 6 => 0, 7 => 1] => [2,1,1,0,1,1,0,1] => [12347]
// array where index is number and value is count of this number in input
void simplify(const int *input, int pocet_pak, int base, int *output, int *count) {
    int *array = (int *) calloc(sizeof(int), pocet_pak); // [0000000]

    for (int i = 0; i < base; i++) {
        int number_of_pak = input[i];
        array[number_of_pak]++;
    }

    *count = 0;
    for (int i = 0; (int) i < pocet_pak; ++i) {
        if (array[i] % 2 == 1) {
            output[*count] = i;
            (*count)++;
        }
    }
}

// [0123456789] 10 => 10 ======= 9 + 80 + 700 + ... + 100000000 = 123456789
int convert_to_number(const int *array, int base, int size) {
    int sum = 0;
    for (int i = 0; i < size; ++i) {
        sum += (int) (array[i] * pow(base, i));
    }
    return sum;
}

// 4 => 01234, 7 => 01234567
//       3210
long get_pak_count(int k) {
    k--;
    long result = 0;
    for (int i = 1; i <= k; ++i) {
        result += (i) * (long) pow(10, k - i);
    }
    return result;
}

int main() {
    clock_t begin = clock();

    // get n and k
    int n, k;
    int get = scanf("%d %d", &n, &k);
    if (get != 2) {
        return 1; //todo: input error
    }
    // get k paks if not 0
    Vector paks[k];
    if (k > 0) {
        for (int i = 0; i < k; ++i) {
            int position_pak;
            get = scanf("%d ", &position_pak);
            if (get != 1) return 0;
            Vector *v = (Vector *) malloc(sizeof(Vector *) + sizeof(int *) * n);
            v->y = position_pak;
            v->vector = (int *) malloc(sizeof(int *) * n);
            for (int j = 0; j < n; j++) {
                char d;
                scanf("%c", &d);
                v->vector[j] = d - 48;
            }
            paks[i] = *v;
        }
    }
    int **maze = init_maze(n);
    int finishX, finishY;
    scanf("%d %d", &finishX, &finishY);
    Point end = {finishX, finishY};
//    printMaze(maze, n);
//    printMaze(apply_pak(paks[0], maze, n), n);
//    printf("%d %d\n", n, k);
//    printf("%d %d\n", finishX, finishY);
//    printMaze(maze, n);
    Result best = check_maze(maze, n, 0, end);
    if (best.size == -1 && k == 0) {
        for (int i = 0; i < n + 2; ++i) {
            free(maze[i]);
        }
        free(maze);
        printf("-1\n");
        clock_t time_end = clock();
        double time_spent = (double) (time_end - begin) / CLOCKS_PER_SEC;
//        printf("%f", time_spent);
        return 0;
    }
    if (best.size == -1) {
        best.size = 10000 * 10000;
    }
    long k_pak = get_pak_count(k) + 1;
    Result *result = (Result *) malloc(sizeof(Result) * k_pak);
    for (int i = 0; i < k_pak; ++i) {
        result[i] = (Result) {-1, -1};
    }

    for (long i = 0; i < (long) pow(k, k) && k != 0; ++i) {
        int pocet_u = 0;
        int *array = convert_base(i, k, &pocet_u); // 21(3) => [2,1,0], 17(3) => [122]
        int *output = (int *) calloc(sizeof(int), k);
        int number = 0;
        simplify(array, k, pocet_u, output, &number);
        if (number == 0) {
            continue;
        }
        int index = convert_to_number(output, k, number);
        int **new_maze = copy_maze(maze, n);
        int max_pak = 0;
        for (int j = 0; j < number; ++j) {
            if (max_pak < paks[output[j]].y) {
                max_pak = paks[output[j]].y;
            }
            new_maze = apply_pak(paks[j], new_maze, n);
        }
        Result r = check_maze(new_maze, n, max_pak, end);
        if (r.size == -1) continue;
        if (result[index].size == -1 || result[index].size > r.size) {
            result[index] = r;
        }

        for (int j = 0; j < pocet_u; ++j) {
            printf("%d", array[j]);
        }
        printf("[%d] =", pocet_u);

        for (int j = 0; j < number; ++j) {
            printf(" %d", output[j]);
        }
        printf("[%d] => %d\n", index, result[index].size);
    }

    for (int i = 0; i < k_pak; ++i) {
        if (result[i].size == -1) continue;
        if (best.size > result[i].size) best = result[i];
    }

    if (best.size != 10000 * 10000) {
        printf("%d\n", best.size);
        for (int i = 0; i < best.count + 1; ++i) {
            if (i == 0) {
                printf("[%d,%d]", k > 0 ? best.way[i]._x : best.way[i]._y, k > 0 ? best.way[i]._y : best.way[i]._x);
            } else {
                printf(",[%d,%d]", k > 0 ? best.way[i]._x : best.way[i]._y, k > 0 ? best.way[i]._y : best.way[i]._x);
            }
        }
    } else {
        printf("-1");
    }


    for (int i = 0; i < n + 2; ++i) {
        free(maze[i]);
    }
    free(maze);
    for (int i = 0; i < k_pak; ++i) {
        free(result[i].way);
    }
    free(result);
    clock_t time_end = clock();
    double time_spent = (double) (time_end - begin) / CLOCKS_PER_SEC;
//    printf("\n%f", time_spent);
    return 0;
}