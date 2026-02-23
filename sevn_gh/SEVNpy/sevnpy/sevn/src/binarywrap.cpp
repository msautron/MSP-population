//
// Created by Giuliano Iorio on 02/09/23.
//

#include "binarywrap.h"
#include "starwrap.h"
#include <binstar.h>
#include <vector>
#include <string>
#include "iowrap.h"
#include "wraputils.h"

thread_local std::vector<size_t> Binstar_property::bpropertiesID_list;
void Binstar_property::fill_bpropertiesID_list(){
    bpropertiesID_list = {
            BWorldtime::ID,
            Semimajor::ID,
            Eccentricity::ID,
            Period::ID,
            RL0::ID,
            RL1::ID,
            GWtime::ID,
            BEvent::ID,
    };
}

//Binary properties implementations
Binstar_property::Binstar_property(npy_intp starting_size) : size{starting_size}, filling_counter{0}{
    //Initialise the static list of binary properties
    if (bpropertiesID_list.empty()) fill_bpropertiesID_list();

    //Initialise all the binary  properties
    for (const auto& propID : bpropertiesID_list){
        bproperties_nparrays[propID]=PyArray_ZEROS(dim,&size,NPY_DOUBLE,0);
    }
    stars_properties.reserve(Binstar_property::nbody);
    //Initialise the two stars
    for (size_t i=0; i<Binstar_property::nbody; i++){
        stars_properties.emplace_back(starting_size);
    }
}
Binstar_property::~Binstar_property(){
    for (const auto& propID : bpropertiesID_list){
        //We are removing this object, IF no one owns the related arrays they can be destroyed
        //so decref them
        Py_DECREF(bproperties_nparrays[propID]);
    }
}
PyObject* Binstar_property::gets(size_t propID, size_t starID){
    return stars_properties[starID].get(propID);
}
void Binstar_property::update_from_binary(Binstar* b){
    //If run out of space, extend the numpy arrays
    if (filling_counter>size-1) extend();
    //Fill the binary arrays
    for  (const auto& propID : bpropertiesID_list){
        insert_value(propID,filling_counter,b->getp(propID));
    }
    //Fill the single stars
    stars_properties[0].update_from_star(b->getstar(0));
    stars_properties[1].update_from_star(b->getstar(1));
    //Update filling counter
    filling_counter++;
}
void Binstar_property::resize(npy_intp new_size){
    PyArray_Dims dims;
    dims.ptr = &new_size;
    dims.len = 1;

    for (auto& prop : bproperties_nparrays){
        auto np_array = prop.second;
        PyObject *NoneObj = PyArray_Resize( (PyArrayObject*)np_array, &dims, 0,  NPY_ANYORDER);
        Py_DECREF(NoneObj);
    }
}
void Binstar_property::extend(){
    size= 2*size+1;
    PyArray_Dims dims;
    dims.ptr = &size;
    dims.len = 1;
    for (auto& prop : bproperties_nparrays){
        auto np_array = prop.second;
        PyObject *NoneObj = PyArray_Resize( (PyArrayObject*)np_array, &dims, 0,  NPY_ANYORDER);
        Py_DECREF(NoneObj);
    }
}

void Binstar_property::trim(){
    //Trim binary
    resize(filling_counter);
    //Trim stars
    stars_properties[0].trim();
    stars_properties[1].trim();
}

Binstar_property _init_binary(Binstar& b){
    Binstar_property bprop;
    //Store the initial value
    bprop.update_from_binary(&b);
    bprop.trim();
    return bprop;
}
Binstar_property _evolve_binary(Binstar& b, double evolve_to){
    Binstar_property bprop;
    bprop.update_from_binary(&b);

    if (evolve_to<0) evolve_to=1E30;
    while(b.getp(BWorldtime::ID)<=evolve_to){
        b.evolve();
        bprop.update_from_binary(&b);
        if (b.breaktrigger()) break;
    }

    bprop.trim();

    return bprop;
}