/**
 * @file initial.h
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Initial functions and variables
 *
 */

// NAN                                                                          
int CHECK_NAN;                                                                  
                                                                                
/* Define pointers */                                                           
double *U, *U0, *U1, *U2, *U3;                                                  
double *Q, *Q0, *Q1, *Q2, *Q3;                                                  
                                                                                
double *U1p, *U1m;                                                              
double *U2p, *U2m;                                                              
                                                                                
double K;                                                                       
                                                                                
/* Define number file */                                                        
int itprint;                                                                    
                                                                                
/* Define freq. output dt and time */                                           
double dtprint, tprint;  

void Allocate_Array();

void New_Size();

void Initial();
