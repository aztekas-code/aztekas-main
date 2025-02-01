/**
 * @file main.h
 * @brief This file contains the Example class declaration.
 */
#ifndef INCLUDE_AZTEKAS_MAIN_H_
#define INCLUDE_AZTEKAS_MAIN_H_

#include <string>

/**
 * Clase que representa un ejemplo básico.
 */
class Example {
public:
  explicit Example(int initial_value);

  void PerformOperation() const;

private:
  int value_; // Valor interno.
};

#endif // INCLUDE_AZTEKAS_MAIN_H_
