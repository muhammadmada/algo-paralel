#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/wait.h>
#include <sys/mman.h>
#include <stdatomic.h>
#include <semaphore.h>

#define MAX_STEPS 100

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
        //DEBUG
        //printf("i: %d, x0: %lf, x1: %lf, y0: %lf, y1: %lf, result: %lf\n", i, x0, x1, y0, y1, result);
        //END DEBUG
    }
    return result;
}

void evaluate(double a, double b, int dx, double *result, double *execution_time) {
    clock_t start_time = clock();
    *result = trapezoidal_integral(a, b, dx);
    clock_t end_time = clock();
    *execution_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    //DEBUG
    //printf("Result: %lf, Execution Time: %lfsec\n", *result, *execution_time);
    //DEBUG END
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

void serial_integration(double a, double b, int dx, double steps, double *result, double *execution_time) {
    printf("--- SERIAL INTEGRATION ---\n");
    for (int step = 1; step <= steps; step++) {
        int n = pow(10, step);
        evaluate(a, b, n, &result[step - 1], &execution_time[step - 1]);

        printf("Step number %d\n", step);
        printf("Result: %lf\n", result[step - 1]);
        printf("Execution time: %lfsec\n\n", execution_time[step - 1]);
    }
    save_to_csv("serial_results.csv", result, execution_time, steps);
    printf("--- SERIAL INTEGRATION END ---\n\n\n");
}

void parallel_integration(double a, double b, int dx, int steps, double *result, double *execution_time) {
    double *shared_results = mmap(NULL, steps * sizeof(double), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    double *shared_exec_time = mmap(NULL, steps * sizeof(double), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

    sem_t *sem = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    sem_init(sem, 1, 0);

    printf("--- PARALLEL INTEGRATION ---\n");

    for (int step = 1; step <= steps; step++) {
        int n = pow(10, step);

        pid_t pid = fork();

        if (pid == 0) {
            double local_result, local_execution_time;
            evaluate(a, n / 2, n, &local_result, &local_execution_time);

            atomic_store_explicit(&shared_results[step - 1], local_result, memory_order_relaxed);
            atomic_store_explicit(&shared_exec_time[step - 1], local_execution_time, memory_order_relaxed);

            sem_post(sem);

            exit(0);
        } else if (pid > 0) {
            double local_result, local_execution_time;
            evaluate((n / 2) + 1, n, n, &local_result, &local_execution_time);

            atomic_store_explicit(&shared_results[step - 1], local_result, memory_order_relaxed);
            atomic_store_explicit(&shared_exec_time[step - 1], local_execution_time, memory_order_relaxed);
        } else {
            perror("Fork error!");
            exit(EXIT_FAILURE);
        }
    }

    for (int step = 1; step <= steps; step++){
        sem_wait(sem);
    }

    for (int step = 0; step < steps; step++) {
        result[step] = atomic_load_explicit(&shared_results[step], memory_order_relaxed);
        execution_time[step] = atomic_load_explicit(&shared_exec_time[step], memory_order_relaxed);
        printf("Step number: %d\nResult: %f\nExecution Time: %fsec\n\n", step + 1, result[step], execution_time[step]);
    }

    printf("--- PARALLEL INTEGRATION END ---\n");

    save_to_csv("parallel_results.csv", result, execution_time, steps);

    munmap(shared_results, steps * sizeof(double));
    munmap(shared_exec_time, steps * sizeof(double));
}



int main() {
    double a = 1.0;
    double b = 5.0;
    int dx = 1000000000;
    double steps = log10(dx);

    double serial_result, serial_execution_time;
    double parallel_result[MAX_STEPS], parallel_execution_time[MAX_STEPS];

    serial_integration(a, b, dx, steps, &serial_result, &serial_execution_time);

    //parallel_integration(a, b, dx, steps, parallel_result, parallel_execution_time);

    return 0;
}