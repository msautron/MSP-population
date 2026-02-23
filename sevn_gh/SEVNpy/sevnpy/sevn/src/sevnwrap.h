/***
 * Created by Giuliano Iorio on 02/09/23.
 * This file contains the proper SEVN wrapper.
 * It contains all the functions and parameters that are then exposed to the Python module +
 * all the necessary to produce the module.
 */
#ifndef SEVN_SEVNWRAP_H
#define SEVN_SEVNWRAP_H

#include <Python.h>
#include <sevn.h>

///SEVNIO HANDLING
static PyObject *sevnio_initialise(PyObject *self, PyObject *args);
static PyObject *sevnio_finalise(PyObject *self, PyObject *args);
static PyObject *sevnio_param(PyObject *self, PyObject *args);
///STELLAR EVOLUTION
static PyObject *evolve_star(PyObject *self, PyObject *args, PyObject *kwargs);
///BINARY EVOLUTION
static PyObject *evolve_binary(PyObject *self, PyObject *args, PyObject *kwargs);

#endif //SEVN_SEVNWRAP_H
