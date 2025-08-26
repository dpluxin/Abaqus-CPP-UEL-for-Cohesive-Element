// ToolKit.h
#pragma once

// ******** Definition of ToolKit funcitons ******** //

// Define the auxilliary function for matrix operations
// Cross product function
static void cross_product(double const *a, double const *b, double *result) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}
// Vector magnitude function
static double vector_magnitude(double const *vec) {
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}
// Normalize vector function
static void normalize(double const *vec, double *result) {
    double mag = vector_magnitude(vec);
    result[0] = vec[0] / mag;
    result[1] = vec[1] / mag;
    result[2] = vec[2] / mag;
}
// End the auxilliary function for matrix operations
