/**
 * @file boundaries.h
 * 
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Boundary condition functions definitions.
 *
 */

void Boundaries(double *B);

void Outflow(double *B);
void Periodic(double *B);
void Reflection(double *B);
void User_Boundaries(double *B);
