/***
 * Created by Giuliano Iorio on 02/09/23.
 * This file contains  utility functions used in the SEVN wrapper
 */

#ifndef SEVN_WRAPUTILS_H
#define SEVN_WRAPUTILS_H

#include <Python.h>
#include <string>
#include <sstream>

std::string PyObject_to_string(PyObject *input);
double PyObject_to_double(PyObject* input);
std::string compose_mass_input(const double& Mass, const std::string& star_flag);


/**
 * Transform a number to a string including the precision.
 * This is a replacement of the to_string method from the standard library
 * that has a constant precision of 6 (so for example it returns 0 for  numbers <1e-7).
 * @param T number
 * @param n precision
 */
template <typename T>
std::string number_to_string_with_precision(const T number, const int n = 15)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << number;
    return std::move(out).str();
}

#endif //SEVN_WRAPUTILS_H
