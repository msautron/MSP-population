//
// Created by Giuliano Iorio on 02/09/23.
//

#include "starwrap.h"
#include <Python.h>
#include <star.h>
#include <vector>
#include <string>
#include "iowrap.h"
#include "wraputils.h"
/***
 * Created by Giuliano Iorio on 02/09/23.
 * This file contains the class Star_property that is used in the wrapper to store the properties
 * of SEVN stars and functions to access to the SEVN single stellar evolution evolution
 */


thread_local std::vector<size_t> Star_property::propertiesID_list;
void Star_property::fill_propertiesID_list(){
    propertiesID_list = {
            Worldtime::ID,
            Localtime::ID,
            Mass::ID,
            Radius::ID,
            Inertia::ID,
            Luminosity::ID,
            Temperature::ID,
            MHE::ID,
            MCO::ID,
            RHE::ID,
            RCO::ID,
            Phase::ID,
            RemnantType::ID,
            PhaseBSE::ID,
            Spin::ID,
            Ebind::ID,
            Zams::ID,
            Event::ID,
            dMRLOdt::ID,
            dMaccwinddt::ID,
            Plife::ID
    };
}

//Star properties implementations
Star_property::Star_property(npy_intp starting_size) : size{starting_size}, filling_counter{0}{
    //Initialise the static list of binary properties
    if (propertiesID_list.empty()) fill_propertiesID_list();

    //Initialise all the properties
    for (const auto& propID : propertiesID_list){
        properties_nparrays[propID]=PyArray_ZEROS(dim,&size,NPY_DOUBLE,0);
    }
}
void Star_property::update_from_star(Star *s){
    //If run out of space, extend the numpy arrays
    if (filling_counter>size-1) extend();
    //Fill the arrays
    for  (const auto& propID : propertiesID_list){
        insert_value(propID,filling_counter,s->getp(propID));
    }
    //Update filling counter
    filling_counter++;
}
void Star_property::resize(npy_intp new_size){
    PyArray_Dims dims;
    dims.ptr = &new_size;
    dims.len = 1;

    for (auto& prop : properties_nparrays){
        auto np_array = prop.second;
        PyObject *NoneObj = PyArray_Resize( (PyArrayObject*)np_array, &dims, 0,  NPY_ANYORDER);
        Py_DECREF(NoneObj);
    }
}
void Star_property::extend(){
    size= 2*size+1;
    PyArray_Dims dims;
    dims.ptr = &size;
    dims.len = 1;
    for (auto& prop : properties_nparrays){
        auto np_array = prop.second;
        PyObject *NoneObj = PyArray_Resize( (PyArrayObject*)np_array, &dims, 0,  NPY_ANYORDER);
        Py_DECREF(NoneObj);
    }
}

Star_property _init_star(Star& s){
    Star_property sprop;
    //Store the initial value
    sprop.update_from_star(&s);
    sprop.trim();
    return sprop;
}
Star_property _evolve_star(Star& s, double evolve_to){

    Star_property sprop;

    //Store the initial value
    sprop.update_from_star(&s);

    if (evolve_to<0) evolve_to=1E30;
    while(s.getp(Worldtime::ID)<=evolve_to){
        s.evolve();
        //if (nstep>sprop.size-1) sprop.extend();
        sprop.update_from_star(&s);
        //nstep++;
        if (s.breaktrigger())  break;
    }

    sprop.trim();

    return sprop;
}
