/**
 * @file main.cpp
 * @brief This file contains the main function of the program.
 */
#include "aztekas/main.h"

/**
 * @brief Main function of the program.
 *
 * Creates an instance of the Example class and performs an operation.
 *
 * @return 0 if the program executes successfully.
 */
int main() {
  Example example(42);

  example.PerformOperation();

  return 0;
}
