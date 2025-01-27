/**
 * @file main.cpp
 * @brief This file contains the main function and the Example class implementation.
 */
#include <iostream>
#include "aztekas/main.h"

Example::Example(int initial_value) : value_(initial_value) {}

void Example::PerformOperation() const {
  std::cout << "Value is: " << value_ << "\n";
}
