//
// Created by Giuliano Iorio on 02/09/23.
//
/***
 * Created by Giuliano Iorio on 02/09/23.
 * This file contains  utility functions used in the SEVN wrapper
 */

#include "wraputils.h"
#include <Python.h>
#include <string>
#include <stdexcept>

std::string PyObject_to_string(PyObject *input){

    //Input owenership borrowed, don't need to DECREF
    std::string output;
    //Check float numbers
    if (PyFloat_Check(input)){
        // Here we use the Python representation that
        // is the string representing  the object in Python
        // so it is already in a string form and we can directly use it
        PyObject* repr_number = PyObject_Repr(input);
        output =  PyUnicode_AsUTF8(repr_number);
        Py_XDECREF(repr_number); // remove the pointer to repr_number
    }
        //Now check Integers or Bool. Notice in Python bool are subclasses of Long, so to disambiguate:
        //Check exact Long (so no Bool)
    else if (PyLong_CheckExact(input)){
        output = std::to_string(PyLong_AsLong(input));
    }
        //Now check Bool True
    else if (PyBool_Check(input) and input==Py_True){
        output = std::string{"true"};
    }
        //Now check Bool False
    else if (PyBool_Check(input) and input==Py_False){
        output = std::string{"false"};
    }
        //Now check the rest of Long subclasses
    else if (PyLong_Check(input)){
        output = std::to_string(PyLong_AsLong(input));
    }
        //Now check strings
    else if (PyUnicode_Check(input)){
        output =  PyUnicode_AsUTF8(input);
    }
        //Oh no, what is this object?
    else{
        throw  std::runtime_error("invalid input in PyObject_to_string");
    }

    return output;
}
double PyObject_to_double(PyObject* input){
    double output;
    if (PyFloat_Check(input))        output = PyFloat_AsDouble(input);
    else if (PyLong_Check(input))    output = PyLong_AsDouble(input);
    else if (PyUnicode_Check(input)){
        auto utf8_str = PyUnicode_AsUTF8(input);
        if (!utf8_str) {
            throw  std::runtime_error("Invalid Unicode input in PyObject_to_double");
        }
        //Convert
        char* endptr;
        output = strtod(utf8_str, &endptr);
        //Check error
        if (endptr == utf8_str) {
            throw  std::runtime_error( "Invalid Unicode double in PyObject_to_double");
        }
    }
    else{
        throw  std::runtime_error("invalid input in PyObject_to_double");
    }

    return output;
}


std::string compose_mass_input(const double& Mass, const std::string& star_flag){

    std::string output;
    output = std::to_string(Mass);

    if (star_flag.empty() or star_flag=="H") return  output;
    else if (star_flag=="HE" or star_flag=="BH" or star_flag=="NS"  or star_flag=="NSEC"
             or star_flag=="HEWD" or star_flag=="COWD" or star_flag=="ONEWD")
        return output + star_flag;
    else{
        std::string err_mess = "Input star_flag "+star_flag+" is invalid\n";
        err_mess+="Valid values are: empty,H,HE,BH,NS,NSEC,HEWD,COWD,ONEWD";
        throw  std::runtime_error(err_mess);
    }
}