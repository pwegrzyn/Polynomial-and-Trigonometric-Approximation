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

double coeffs[NUMBER_OF_BASE_FUNCTIONS];
double approx_nodes_x[NUMBER_OF_NODES];
double approx_nodes_y[NUMBER_OF_NODES];
double intermediate_points_x[NUMBER_OF_NODES - 1];
double intermediate_points_y[NUMBER_OF_NODES - 1];
double approximated_points_y[NUMBER_OF_NODES - 1];
double approximation_errors[NUMBER_OF_NODES - 1];

void gaussian_elimination(double *a, double *b, double *x, int n);
void gen_approx_nodes(double (*f)(double));
void polynomial_approximation(double (*f)(double));
double f1(double x);
double f2(double x);
double eval_approx_poly(double x);
void print_results_poly();
void gen_plot_data_poly(FILE *file, FILE *file_int, double (*f)(double));
void approximate_points_poly(double (*f)(double));

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

void polynomial_approximation(double (*f)(double)) {
    double *g = (double*)malloc(NUMBER_OF_BASE_FUNCTIONS * NUMBER_OF_BASE_FUNCTIONS * sizeof(double));
    double *b = (double*)malloc(NUMBER_OF_BASE_FUNCTIONS * sizeof(double));
    if(!g || !b) {
        perror("Error while allocating memory\n");
    }
    // Generate G matrix
    for(int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++) {
        for(int j = 0; j < NUMBER_OF_BASE_FUNCTIONS; j++) {
            for(int k = 0; k < NUMBER_OF_NODES; k++) {
                g[i * NUMBER_OF_BASE_FUNCTIONS + j] += pow(approx_nodes_x[k], i + j);
            }
        }
    }
    // Generate B matrix
    for(int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++) {
        for(int j = 0; j < NUMBER_OF_NODES; j++) {
            b[i] += f(approx_nodes_x[j]) * pow(approx_nodes_x[j], i);
        }
        coeffs[i] = 0;
    }
    // G * coeffs = B
    gaussian_elimination(g, b, coeffs, NUMBER_OF_BASE_FUNCTIONS);
    free(g);
    free(b);
}

// [-2pi, 2pi]
double f1(double x) {
    return 10 + pow(x, 2) * 0.5 - 10 * cos(2 * x);
}

// [-100, 100]
double f2(double x) {
    return -x * sin(sqrt(3 * fabs(x - 1)));
}

void approximate_points_poly(double (*f)(double)) {
    gen_approx_nodes(f);
    polynomial_approximation(f);
    for(int i = 0; i < NUMBER_OF_NODES - 1; i++) {
        approximated_points_y[i] = eval_approx_poly(intermediate_points_x[i]);
        approximation_errors[i] = fabs(approximated_points_y[i] - intermediate_points_y[i]) / fabs(intermediate_points_y[i]);
    }
}

double eval_approx_poly(double x) {
    double res = 0;
    for(int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++) {
        res += coeffs[i] * pow(x, i);
    }
    return res;
}

void print_results_poly() {
    printf("Apporximation nodes:\n");
    for(int i = 0; i < NUMBER_OF_NODES; i++) {
        printf("    x = %10f  |  y = %10f\n", approx_nodes_x[i], approx_nodes_y[i]);
    }
    printf("Approximation polynomial coefficients:\n");
    for(int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++) {
        printf("    a%d = %10f |\n", i, coeffs[i]);
    }
    printf("Approximation in intermediate points results:\n");
    for(int i = 0; i < NUMBER_OF_NODES - 1; i++) {
        printf("    x = %10f  |   real_y = %10f  |  approx_y = %10f  |  relative_error = %10f\n", 
            intermediate_points_x[i], intermediate_points_y[i], approximated_points_y[i], approximation_errors[i]);
    }
    printf("----------------------------------------------------\n\n");
}

void gen_plot_data_poly(FILE *file, FILE *file_int, double (*f)(double)) {
    double h = ((double)DOMAIN_SIZE) / (PLOT_GRANULITY - 1);
    char buffer[100];
    for(int i = 0; i < PLOT_GRANULITY; i++) {
        double x = ((double)i) * h + DOMAIN_START;
        double y = eval_approx_poly(x);
        sprintf(buffer, "%f, %f\n", x, y);
        fwrite(buffer, sizeof(char), strlen(buffer), file_int);
        
        x = ((double)i) * h + DOMAIN_START;
        y = f(x);
        sprintf(buffer, "%f, %f\n", x, y);
        fwrite(buffer, sizeof(char), strlen(buffer), file);
    }
}

void gaussian_elimination(double *a, double *b, double *x, int n) {

    int column, row, diagonal, max_pivot_row, j;
    double max_pivot, tmp;
    for(diagonal = 0; diagonal < n; diagonal++) {
        max_pivot_row = diagonal;
        max_pivot = *(a + (diagonal * n + diagonal));    // i,ith element of the matrix
        for(row = diagonal + 1; row < n; row++) {
            tmp = fabs(*(a + (row * n + diagonal)));
            if(tmp > max_pivot) {
                max_pivot_row = row;
                max_pivot = tmp;
            }
        }

        if(diagonal != max_pivot_row) {
            for(int k = 0; k < n; k++) {
                double *tmp_pointer1 = a + (diagonal * n + k);
                double *tmp_pointer2 = a + (max_pivot_row * n + k);
                tmp = *tmp_pointer1;
                *tmp_pointer1 = *tmp_pointer2;
                *tmp_pointer2 = tmp;
            }
            tmp = b[diagonal];
            b[diagonal] = b[max_pivot_row];
            b[max_pivot_row] = tmp;
        }

        for(row = diagonal + 1; row < n; row++) {
            tmp = *(a + (row * n + diagonal)) / *(a + (diagonal * n + diagonal));
            for(column = diagonal + 1; column < n; column++) {
                *(a + (row * n + column)) -= tmp * *(a + (diagonal * n + column));
            }
            *(a + (row * n + diagonal)) = 0;
            b[row] -= tmp * b[diagonal];
        }
    }

    for(row = n - 1; row >= 0; row--) {
        tmp = b[row];
        for(j = n - 1; j > row; j--) {
            tmp -= x[j] * *(a + (row * n + j));
        }
        x[row] = tmp / *(a + (row * n + row));
    }

}

int main() {

    FILE *data;
    FILE *approx_data;
    data = fopen("data.csv", "w+");
    approx_data = fopen("approx_data.csv", "w+");

    approximate_points_poly(f2);
    print_results_poly();
    gen_plot_data_poly(data, approx_data, f2);

}