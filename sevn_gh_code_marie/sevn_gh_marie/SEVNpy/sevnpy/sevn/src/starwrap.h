/***
 * Created by Giuliano Iorio on 02/09/23.
 * This file contains the class Star_property that is used in the wrapper to store the properties
 * of SEVN stars and functions to access to the SEVN single stellar evolution evolution
 */

#ifndef SEVN_STARWRAP_H
#define SEVN_STARWRAP_H

/** Include Python **/
#include <Python.h>
/** Include Numpy API **/
/*! @brief If the extension is splitted in multiple files, we have to add the following definitions
 * (https://docs.scipy.org/doc/numpy-1.10.1/reference/c-api.array.html#miscellaneous) */
#define NO_IMPORT_ARRAY // To be add in files using the NUMPY API but not directlly calling import_array()
#define PY_ARRAY_UNIQUE_SYMBOL SEVNWRAP_ARRAY_API
/*! @brief Use the NumPy 1.7 API. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
/** Other includes **/
#include <vector>
#include <string>
#include <map>
#include <iostream>

class Star;


class Star_property{

public:

    Star_property(npy_intp starting_size=600);

    ~Star_property(){
        for (const auto& propID : propertiesID_list){
            //We are removing this object, IF no one owns the related arrays they can be destroyed
            //so decref them
            Py_DECREF(properties_nparrays[propID]);
        }
    }

    //This vector contains the ID of all the properties that will be considered in output
    //It is a static vector that is based on static variables in the sevn static_main.h
    //the order of static variable initialisation coming from different translation units
    //is undefined, so we cannot be sure that the sevn ID have been already initiliased
    //for this reason we initialise an empty vector in the cpp file and
    //we use a method (fill_propertiesID_list) called by the constructor to
    //initialise it at runtime if the vector is empty (so it is initialised just once per thread)
    static thread_local std::vector<size_t> propertiesID_list;

    //Get the pointer to the PyObject storing the given property
    inline PyObject* get(size_t propID){
        return properties_nparrays[propID];
    }

    //Fill the array from with properties from  a star
    void update_from_star(Star* s);

    //Trim the unnecessary space
    inline void trim(){
        resize(filling_counter);
    }

private:

    const npy_intp dim =1; //Array dimension is 1 (1D Array)
    npy_intp size =600; //Default starting length
    npy_intp filling_counter = 0;

    //This vector contains a map from the SEVN property ID to the sequential ID of the list of numpy arrays storing the properties
    //It is needed because it is not assured the ID in the propertiesID_list  are sorted and sequential
    std::map<size_t,PyObject*> properties_nparrays;

    //Insert a new value for a given property
    inline void insert_value(size_t propID, npy_intp index, double value){
        *reinterpret_cast<double *>(PyArray_GETPTR1(reinterpret_cast<PyArrayObject *>(properties_nparrays[propID]), index)) = value;
    }

    void resize(npy_intp new_size);

    void extend();

    void fill_propertiesID_list();

};

Star_property _init_star(Star& s);
Star_property _evolve_star(Star& s, double evolve_to = -1.0);

#endif //SEVN_STARWRAP_H
