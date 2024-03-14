#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#define MAX_STEPS 100

typedef struct {
    double a;
    double b;
    int dx;
    int step;
    double *result;
    double *execution_time;
} ThreadData;

typedef struct {
    pthread_t thread_id;
    int step;
} ThreadInfo;

double function(double x) {
    return x * x;
}

double trapezoidal_integral(double a, double b, int dx) {
    double h = (b - a) / dx;
    double result = 0.0;
    for (int i = 0; i < dx; i++) {
        double x0 = a + i * h;
        double x1 = x0 + h;

        double y0 = function(x0);
        double y1 = function(x1);

        result += 0.5 * (y0 + y1) * h;
    }
    return result;
}

void *evaluate_step(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    clock_t start_time = clock();
    data->result[data->step] = trapezoidal_integral(data->a, data->b, data->dx);
    clock_t end_time = clock();
    data->execution_time[data->step] = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    return NULL;
}

void save_to_csv(const char *filename, double *result, double *execution_time, int steps) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "Step,Result,Execution Time\n");
    for (int step = 0; step < steps; step++) {
        fprintf(file, "%d,%lf,%lf\n", step + 1, result[step], execution_time[step]);
    }

    fclose(file);
}

void parallel_integration(double a, double b, int dx, int steps, double *result, double *execution_time, int *num_threads) {
    ThreadData thread_data[MAX_STEPS];
    ThreadInfo thread_info[MAX_STEPS];

    printf("--- PARALLEL INTEGRATION ---\n");

    *num_threads = 0;

    for (int step = 1; step <= steps; step++) {
        thread_data[step - 1].a = a;
        thread_data[step - 1].b = b;
        thread_data[step - 1].dx = pow(10, step);
        thread_data[step - 1].step = step - 1;
        thread_data[step - 1].result = result;
        thread_data[step - 1].execution_time = execution_time;

        if (pthread_create(&thread_info[step - 1].thread_id, NULL, evaluate_step, &thread_data[step - 1]) != 0) {
            perror("Error creating thread");
            exit(EXIT_FAILURE);
        }

        thread_info[step - 1].step = step;
        (*num_threads)++;
    }

    for (int step = 1; step <= steps; step++) {
        if (pthread_join(thread_info[step - 1].thread_id, NULL) != 0) {
            perror("Error joining thread");
            exit(EXIT_FAILURE);
        }

        printf("Thread %d processed step %d\n", (int)thread_info[step - 1].thread_id, thread_info[step - 1].step);
        printf("Result: %f\nExecution Time: %fsec\n\n", result[step - 1], execution_time[step - 1]);
    }

    printf("--- PARALLEL INTEGRATION END ---\n");

    save_to_csv("parallel_results.csv", result, execution_time, steps);
}

void serial_integration(double a, double b, int dx, double steps, double *result, double *execution_time) {
    printf("--- SERIAL INTEGRATION ---\n");
    for (int step = 1; step <= steps; step++) {
        int n = pow(10, step);
        ThreadData data = {a, b, n, step - 1, result, execution_time};
        evaluate_step(&data);
        printf("Step number %d\n", step);
        printf("Result: %lf\n", result[step - 1]);
        printf("Execution time: %lfsec\n\n", execution_time[step - 1]);
    }
    save_to_csv("serial_results.csv", result, execution_time, steps);
    printf("--- SERIAL INTEGRATION END ---\n\n\n");
}

int main() {
    double a = 1.0;
    double b = 5.0;
    int dx = 1000000000;
    double steps = log10(dx);

    double parallel_result[MAX_STEPS], parallel_execution_time[MAX_STEPS];
    double serial_result[MAX_STEPS], serial_execution_time[MAX_STEPS];
    int num_threads;

    serial_integration(a, b, dx, steps, serial_result, serial_execution_time);
    parallel_integration(a, b, dx, steps, parallel_result, parallel_execution_time, &num_threads);

    printf("Number of threads used in parallel integration: %d\n", num_threads);

    return 0;
}
