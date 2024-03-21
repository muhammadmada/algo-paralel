#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

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

void parallel_integration(double a, double b, int dx, long int steps, double *result, double *execution_time) {
    printf("--- PARALLEL INTEGRATION ---\n");

    double local_result[MAX_STEPS] = {0};
    double local_execution_time[MAX_STEPS] = {0};

    omp_set_num_threads(6);

    #pragma omp parallel for
    for (long int step = 1; step <= steps; step++) {
        
        double step_result, step_execution_time;
        evaluate(a, b, dx, &step_result, &step_execution_time);

        
        local_result[step - 1] = step_result;
        local_execution_time[step - 1] = step_execution_time;
    }

    
    #pragma omp parallel for
    for (long int i = 0; i < steps; i++) {
        #pragma omp atomic
        result[i] += local_result[i];
        #pragma omp atomic
        execution_time[i] += local_execution_time[i];
        printf("Step number %li\n", i);
        printf("Result: %lf\n", result[i - 1]);
        printf("Execution time: %lf sec\n\n", execution_time[i - 1]);
    }

    printf("--- PARALLEL INTEGRATION END ---\n");

    save_to_csv("parallel_results.csv", result, execution_time, steps);
}

void serial_integration(double a, double b, int dx, long int steps, double *result, double *execution_time) {
    printf("--- SERIAL INTEGRATION ---\n");
    for (long int step = 1; step <= steps; step++) {
        int n = pow(10, step);
        evaluate(a, b, dx, &result[step - 1], &execution_time[step - 1]);
        printf("Step number %li\n", step);
        printf("Result: %lf\n", result[step - 1]);
        printf("Execution time: %lf sec\n\n", execution_time[step - 1]);
    }
    save_to_csv("serial_results.csv", result, execution_time, steps);
    printf("--- SERIAL INTEGRATION END ---\n\n\n");
}

int main() {
    double a = 1.0;
    double b = 5.0;
    int dx = 1000000000;
    long int steps = log10(dx);

    double parallel_result[MAX_STEPS] = {0}, parallel_execution_time[MAX_STEPS] = {0};
    double serial_result[MAX_STEPS], serial_execution_time[MAX_STEPS];
    int num_threads;

    //serial_integration(a, b, dx, steps, serial_result, serial_execution_time);
    parallel_integration(a, b, dx, steps, parallel_result, parallel_execution_time);

    return 0;
}
