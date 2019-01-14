#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define PI 3.14157
#define NUMBER_OF_NODES 20
#define NUMBER_OF_BASE_FUNCTIONS 20
#define PLOT_GRANULITY 1000
#define DOMAIN_SIZE 200
#define DOMAIN_START -100

double coeffs_a[NUMBER_OF_BASE_FUNCTIONS];
double coeffs_b[NUMBER_OF_BASE_FUNCTIONS];
double approx_nodes_x[NUMBER_OF_NODES];
double approx_nodes_y[NUMBER_OF_NODES];
double intermediate_points_x[NUMBER_OF_NODES - 1];
double intermediate_points_y[NUMBER_OF_NODES - 1];
double approximated_points_y[NUMBER_OF_NODES - 1];
double approximation_errors[NUMBER_OF_NODES - 1];

void gen_approx_nodes(double (*f)(double));
void trigonometric_approximation(double (*f)(double));
double f1(double x);
double f2(double x);
double eval_approx_trig(double x);
void print_results_trig();
void gen_plot_data_trig(FILE *file, FILE *file_int, double (*f)(double));
void approximate_points_trig(double (*f)(double));

void gen_approx_nodes(double (*f)(double)) {
    double h = ((double)DOMAIN_SIZE) / (NUMBER_OF_NODES - 1);
    for(int i = 0; i < NUMBER_OF_NODES; i++) {
        approx_nodes_x[i] = ((double)i) * h + DOMAIN_START;
        approx_nodes_y[i] = f(approx_nodes_x[i]);
    }
    double start =  h / 2; 
    for(int i = 0; i < NUMBER_OF_NODES - 1; i++) {
        intermediate_points_x[i] = start + ((double)i) * h + DOMAIN_START;
        intermediate_points_y[i] = f(intermediate_points_x[i]);
    }
}

void trigonometric_approximation(double (*f)(double)) {
    for(int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++) {
        coeffs_a[i] = 0;
        coeffs_b[i] = 0;
        for(int j = 0; j < NUMBER_OF_NODES; j++) {
            coeffs_a[i] += f(approx_nodes_x[j]) * cos((i * approx_nodes_x[j] * PI) / (DOMAIN_SIZE / 2));
            coeffs_b[i] += f(approx_nodes_x[j]) * sin((i * approx_nodes_x[j] * PI) / (DOMAIN_SIZE / 2));
        }
        coeffs_a[i] /= NUMBER_OF_NODES / 2;
        coeffs_b[i] /= NUMBER_OF_NODES / 2;
    }
}

// [-2pi, 2pi]
double f1(double x) {
    return 10 + pow(x, 2) * 0.5 - 10 * cos(2 * x);
}

// [-100, 100]
double f2(double x) {
    return -x * sin(sqrt(3 * fabs(x - 1)));
}

void approximate_points_trig(double (*f)(double)) {
    gen_approx_nodes(f);
    trigonometric_approximation(f);
    for(int i = 0; i < NUMBER_OF_NODES - 1; i++) {
        approximated_points_y[i] = eval_approx_trig(intermediate_points_x[i]);
        approximation_errors[i] = fabs(approximated_points_y[i] - intermediate_points_y[i]) / fabs(intermediate_points_y[i]);
    }
}

double eval_approx_trig(double x) {
    double res = 0;
    for(int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++) {
        if(i == 0) {
            res += 0.5 * coeffs_a[i];
        } else {
            res += coeffs_a[i] * cos((i * x * PI) / (DOMAIN_SIZE / 2)) + coeffs_b[i] * sin((i * x * PI) / (DOMAIN_SIZE / 2));
        }
    }
    return res;
}

void print_results_trig() {
    printf("Apporximation nodes:\n");
    for(int i = 0; i < NUMBER_OF_NODES; i++) {
        printf("    x = %10f  |  y = %10f\n", approx_nodes_x[i], approx_nodes_y[i]);
    }
    printf("Approximation trig functions coefficients:\n");
    for(int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++) {
        printf("    a%d = %10f |  b%d = %10f\n", i, coeffs_a[i], i, coeffs_b[i]);
    }
    printf("Approximation in intermediate points results:\n");
    for(int i = 0; i < NUMBER_OF_NODES - 1; i++) {
        printf("    x = %10f  |   real_y = %10f  |  approx_y = %10f  |  relative_error = %10f\n", 
            intermediate_points_x[i], intermediate_points_y[i], approximated_points_y[i], approximation_errors[i]);
    }
    printf("----------------------------------------------------\n\n");
}

void gen_plot_data_trig(FILE *file, FILE *file_int, double (*f)(double)) {
    double h = ((double)DOMAIN_SIZE) / (PLOT_GRANULITY - 1);
    char buffer[100];
    for(int i = 0; i < PLOT_GRANULITY; i++) {
        double x = ((double)i) * h + DOMAIN_START;
        double y = eval_approx_trig(x);
        sprintf(buffer, "%f, %f\n", x, y);
        fwrite(buffer, sizeof(char), strlen(buffer), file_int);
        
        x = ((double)i) * h + DOMAIN_START;
        y = f(x);
        sprintf(buffer, "%f, %f\n", x, y);
        fwrite(buffer, sizeof(char), strlen(buffer), file);
    }
}

int main() {

    FILE *data;
    FILE *approx_data;
    data = fopen("data.csv", "w+");
    approx_data = fopen("approx_data.csv", "w+");

    approximate_points_trig(f2);
    print_results_trig();
    gen_plot_data_trig(data, approx_data, f2);

}