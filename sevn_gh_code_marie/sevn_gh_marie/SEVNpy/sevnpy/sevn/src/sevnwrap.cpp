/***
 * Created by Giuliano Iorio on 02/09/23.
 * This file contains the proper SEVN wrapper.
 * It contains all the functions and parameters that are then exposed to the Python module +
 * all the necessary to produce the module.
 */

/** Other Python **/
#include <Python.h>
/** Include Numpy API **/
/*! @brief If the extension is splitted in multiple files, we have to add the following definitions
 * (https://docs.scipy.org/doc/numpy-1.10.1/reference/c-api.array.html#miscellaneous) */
#define PY_ARRAY_UNIQUE_SYMBOL SEVNWRAP_ARRAY_API
/*! @brief Use the NumPy 1.7 API. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
/** Other includes **/
#include <vector>
#include <string>
/** Wrapper includes **/
#include "sevnwrap.h"
#include "starwrap.h"
#include "binarywrap.h"
#include "iowrap.h"
#include "wraputils.h"



///SEVNIO HANDLING
static PyObject *sevnio_finalise(PyObject *self, PyObject *args) {
    //First of all check if the IO_handler is initialised
    IO_handler* sevnio = check_and_return_IOinstance();
    IO_handler::close();
    Py_RETURN_NONE;
}
static PyObject *sevnio_initialise(PyObject *self, PyObject *args){

    try{
        PyObject* dict = nullptr;
        PyArg_ParseTuple(args, "|O", &dict);
        IO_handler::initialise(dict);
    }
        //Catch possible initialisation errors and throw a Python error
    catch (const std::exception& err){ //Catch all errors
        PyErr_SetString(PyExc_RuntimeError,err.what());
        Py_RETURN_NONE;
    }

    Py_RETURN_NONE;
}
static PyObject *sevnio_param(PyObject *self, PyObject *args){

    IO_handler* iohandl;
    IO* sevnio;
    SEVNpar* svpar = nullptr;
    bool to_be_deleted = false;


    //Try to get the instance, if not initialised will return an error,
    //so catch it and set again to nullptr
    try{
        iohandl = IO_handler::getInstance();
        sevnio=iohandl->get_sevnio();
        svpar = &sevnio->svpar;
    } catch (std::runtime_error& err){
        iohandl = nullptr;
        sevnio = nullptr;
    };

    //If not loaded
    if (svpar == nullptr){
        sevnio = new IO();
        svpar = &sevnio->svpar;
        to_be_deleted=true;
    }
    svpar->set_myself();




    PyObject* outputDictionary = PyDict_New();
    PyObject* auxiliaryDictionary = PyDict_New();
    //Anyway get the bool map
    auto bmap = svpar->get_bool_map();
    for (auto v: bmap){
        const auto& name        = v.first.c_str();
        const bool& val         = v.second.first;
        const auto& description = v.second.second.c_str();
        PyDict_SetItemString(outputDictionary,name,PyBool_FromLong(val));
        PyDict_SetItemString(auxiliaryDictionary,name,PyUnicode_FromString(description));
    }

    auto smap = svpar->get_str_map();
    for (auto v: smap){
        const auto& name        = v.first.c_str();
        const auto& val         = v.second.first.c_str();
        const auto& description = v.second.second.c_str();
        PyDict_SetItemString(outputDictionary,name,PyUnicode_FromString(val));
        PyDict_SetItemString(auxiliaryDictionary,name,PyUnicode_FromString(description));
    }

    auto nmap = svpar->get_num_map();
    for (auto v: nmap){
        const auto& name          = v.first.c_str();
        const double& val         = v.second.first;
        const auto& description   = v.second.second.c_str();
        PyDict_SetItemString(outputDictionary,name,PyFloat_FromDouble(val));
        PyDict_SetItemString(auxiliaryDictionary,name,PyUnicode_FromString(description));
    }


    if (to_be_deleted){
        delete sevnio;
        sevnio = nullptr;
    }

    return PyTuple_Pack(2,outputDictionary,auxiliaryDictionary);
}

///STELLAR EVOLUTION
static PyObject *evolve_star(PyObject *self, PyObject *args, PyObject *kwargs){

    IO_handler* sevnio;
    //First of all check if the IO_handler is initialised
    try{
        sevnio = IO_handler::getInstance();
        if (sevnio== nullptr) throw std::runtime_error("SEVN manager has not been initiliased, please run sevnio_initialise() before calling evolve_star");
    }
    catch (const std::exception& err){ //Catch all errors
        PyErr_SetString(PyExc_RuntimeError,err.what());
        Py_RETURN_NONE;
    }


    //Input parameters
    static char *kwlist[] = {strdup("Mzams"),
                             strdup("Z"),
                             strdup("spin"),
                             strdup("tstart"),
                             strdup("tend"),
                             strdup("star_flag"),
                             strdup("snmodel"),
                             strdup("rseed"),
                             strdup("Mass"),
                             strdup("MHE"),
                             strdup("MCO"),
                             strdup("Radius"),
                             strdup("just_init"),
                             NULL
    };

    //Define the parameter to be read
    double Mzams,Z;
    double spin{0};
    bool just_init{false};
    unsigned long rseed{0};
    const char*  star_flagCstring = "H";
    const char*  snmodelCstring = "rapid";



    //Define the special input (they can be numbers or strings)
    PyObject* tstart = NULL; //PyUnicode_FromString("zams");
    PyObject* tend   = NULL; //PyUnicode_FromString("end");

    PyObject* Mcurrent          = NULL;
    PyObject* MHEcurrent        = NULL;
    PyObject* MCOcurrent        = NULL;
    PyObject* Rcurrent          = NULL;

    //Parse
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dd|d$OOsskOOOOp",kwlist,
                                     &Mzams,
                                     &Z,
                                     &spin,
                                     &tstart,
                                     &tend,
                                     &star_flagCstring,
                                     &snmodelCstring,
                                     &rseed,
                                     &Mcurrent,
                                     &MHEcurrent,
                                     &MCOcurrent,
                                     &Rcurrent,
                                     &just_init)){
        return nullptr;
    }
    //NOTICE: PyArg_ParseTupleAndKeywords does not increase the reference counter of the PyObject
    //So we don't have to DECREF tstart and tend.
    //Notice also that if we check the ref value of tstart and tend here can be much higher than
    //2 because it reflects the reference counter in the main python script.




    /**** DEFINE THE INPUT VECTOR *****/
    std::string _mzams,_z,_spin;
    std::string _tstart{"zams"},_tend{"end"};
    std::string star_flag(star_flagCstring),_snmodel(snmodelCstring);

    //Compose the Mzams input (it can have suffixes) taking into account errors
    try{
        _mzams   = compose_mass_input(Mzams,star_flag);
    }
    catch (const std::exception& err){ //Catch all errors
        PyErr_SetString(PyExc_RuntimeError,err.what());
        Py_RETURN_NONE;
    }

    _z       = number_to_string_with_precision(Z,15);
    _spin    = number_to_string_with_precision(spin,15);

    //Convert Pyobject to string taking into caoount errors
    try{
        if(tstart) _tstart = PyObject_to_string(tstart);
        if(tend)   _tend   = PyObject_to_string(tend);
    } catch (const std::exception& err){ //Catch all errors
        PyErr_SetString(PyExc_RuntimeError,err.what());
        Py_RETURN_NONE;
    }


    //Set the input parameters
    std::vector<std::string> init_params{_mzams,_z,_spin,_snmodel,_tstart,_tend,"all"};
    //Set fake IDs
    size_t ID = 0;

    //Initialise star and evolve taking into account errors
    Star_property *sprop;
    std::string star_name;
    long star_rseed=0;
    double star_tlife;
    try{
        //Initialise the star
        Star s1(sevnio->get_sevnio(),init_params,ID,false,rseed);
        star_name  = s1.get_name();
        star_rseed = s1.get_rseed();

        //Update value
        if (Mcurrent and Mcurrent!=Py_None){
            double _Mass = PyObject_to_double(Mcurrent);
            s1.update_from_binary(Mass::ID,_Mass-s1.getp(Mass::ID));
        }
        if (MHEcurrent and MHEcurrent!=Py_None){
            double _MHE = PyObject_to_double(MHEcurrent);
            s1.update_from_binary(MHE::ID,_MHE-s1.getp(MHE::ID));
        }
        if (MCOcurrent and MCOcurrent!=Py_None){
            double _MCO = PyObject_to_double(MCOcurrent);
            s1.update_from_binary(MCO::ID,_MCO-s1.getp(MCO::ID));
        }
        if (Rcurrent and Rcurrent!=Py_None){
            double _RCO = PyObject_to_double(Rcurrent);
            s1.update_from_binary(Radius::ID,_RCO-s1.getp(Radius::ID));
        }

        //If just init just initialise to tini
        if (just_init) sprop = new Star_property(_init_star(s1));
        else sprop = new Star_property(_evolve_star(s1));

        //Estimate tlife post evolution (because maybe meanwhile it changed trakcs)
        star_tlife = s1.get_tphase()[Lookup::Phases::Remnant];
    }
    catch (const std::exception& err){ //Catch all errors
        PyErr_SetString(PyExc_RuntimeError,err.what());
        Py_RETURN_NONE;
    }




    //Create the output Dictionary containing the properties
    PyObject* outputProperyDictionary = PyDict_New();
    for (const auto& propID : Star_property::propertiesID_list){
        PyDict_SetItemString(outputProperyDictionary, Property::all[propID]->name().c_str(), sprop->get(propID));
    }

    //Create the output Dictionary containing the additional informations
    PyObject* outputAuxiliaryDictionary = PyDict_New();
    PyDict_SetItemString(outputAuxiliaryDictionary,"Log",PyUnicode_FromString(sevnio->get_log().c_str()));
    PyDict_SetItemString(outputAuxiliaryDictionary,"rseed",PyLong_FromLong(star_rseed));
    PyDict_SetItemString(outputAuxiliaryDictionary,"ID",PyLong_FromSsize_t(ID));
    PyDict_SetItemString(outputAuxiliaryDictionary,"name",PyUnicode_FromString(star_name.c_str()));
    PyDict_SetItemString(outputAuxiliaryDictionary,"tlife",PyFloat_FromDouble(star_tlife));

    //Reset log
    sevnio->reset_log();
    //Clean sprop
    delete sprop;
    sprop = nullptr;

    PyObject* return_tuple = PyTuple_Pack(2,outputProperyDictionary,outputAuxiliaryDictionary);
    //PyTuple_Pack increments by one the ref count, so just decref it to maintain a value of 1
    Py_DECREF(outputProperyDictionary);
    Py_DECREF(outputAuxiliaryDictionary);

    return return_tuple;

}

///BINSTAR EVOLUTION
static PyObject *evolve_binary(PyObject *self, PyObject *args, PyObject *kwargs){
    IO_handler* sevnio;
    //First of all check if the IO_handler is initialised
    try{
        sevnio = IO_handler::getInstance();
        if (sevnio== nullptr) throw std::runtime_error("SEVN manager has not been initiliased, "
                                                       "please run sevnio_initialise() before calling evolve_binary");
    }
    catch (const std::exception& err){ //Catch all errors
        PyErr_SetString(PyExc_RuntimeError,err.what());
        Py_RETURN_NONE;
    }

    //Input parameters
    static char *kwlist[] = {strdup("Semimajor"),
                             strdup("Eccentricity"),
                             strdup("Mzams_0"),
                             strdup("Z_0"),
                             strdup("Mzams_1"),
                             strdup("Z_1"),
                             strdup("spin_0"),
                             strdup("tstart_0"),
                             strdup("spin_1"),
                             strdup("tstart_1"),
                             strdup("tend"),
                             strdup("snmodel"),
                             strdup("star_flag_0"),
                             strdup("star_flag_1"),
                             strdup("rseed"),
                             strdup("Mass_0"),
                             strdup("MHE_0"),
                             strdup("MCO_0"),
                             strdup("Radius_0"),
                             strdup("Mass_1"),
                             strdup("MHE_1"),
                             strdup("MCO_1"),
                             strdup("Radius_1"),
                             strdup("just_init"),
                             NULL
    };

    //Define the parameter to be read
    double semimajor{0.}, eccentricity{0.};
    double Mzams[2]={0.,0.};
    double Z[2]={0.,0.};
    double spin[2]={0.,0.};
    bool just_init{false};
    unsigned long rseed{0};
    const char*  star_flagCstring[2] = {"H","H"};
    //const char*  star_flagCstring_1 = "H";
    const char*  snmodelCstring = "rapid";

    //Define the special input (they can be numbers or strings)
    PyObject* tstart[2] = {NULL,NULL}; //PyUnicode_FromString("zams");
    PyObject* tend   = NULL; //PyUnicode_FromString("end");

    PyObject* Mcurrent[2]          = {NULL,NULL};
    PyObject* MHEcurrent[2]        = {NULL,NULL};
    PyObject* MCOcurrent[2]        = {NULL,NULL};
    PyObject* Rcurrent[2]          = {NULL,NULL};

    //Parse
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dddddd|$dOdOOssskOOOOOOOOp",kwlist,
                                     &semimajor,
                                     &eccentricity,
                                     &Mzams[0],
                                     &Z[0],
                                     &Mzams[1],
                                     &Z[1],
                                     &spin[0],
                                     &tstart[0],
                                     &spin[1],
                                     &tstart[1],
                                     &tend,
                                     &snmodelCstring,
                                     &star_flagCstring[0],
                                     &star_flagCstring[1],
                                     &rseed,
                                     &Mcurrent[0],
                                     &MHEcurrent[0],
                                     &MCOcurrent[0],
                                     &Rcurrent[0],
                                     &Mcurrent[1],
                                     &MHEcurrent[1],
                                     &MCOcurrent[1],
                                     &Rcurrent[1],
                                     &just_init)){
        return nullptr;
    }
    /**** DEFINE THE INPUT VECTOR *****/
    std::string _semimajor,_eccentricity;
    std::string _mzams[2],_z[2],_spin[2];
    std::string _tstart[2]{"zams","zams"};
    std::string _tend{"end"};
    std::string star_flag[2]={star_flagCstring[0],star_flagCstring[1]};
    std::string _snmodel(snmodelCstring);


    //Global properties
    _semimajor    = number_to_string_with_precision(semimajor,15);
    _eccentricity = number_to_string_with_precision(eccentricity,15);
    try {
        if (tend) _tend = PyObject_to_string(tend);
    } catch (const std::exception& err){ //Catch all errors
        PyErr_SetString(PyExc_RuntimeError,err.what());
        Py_RETURN_NONE;
    }
    //Star properties
    for (size_t i=0; i<2; i++){
        //Compose the Mzams input (it can have suffixes) taking into account errors
        try{
            _mzams[i]   = compose_mass_input(Mzams[i],star_flag[i]);
        }
        catch (const std::exception& err){ //Catch all errors
            PyErr_SetString(PyExc_RuntimeError,err.what());
            Py_RETURN_NONE;
        }
        //Z and spin
        _z[i]       = number_to_string_with_precision(Z[i],15);
        _spin[i]    = number_to_string_with_precision(spin[i],15);
        //Convert Pyobject to string taking into caoount errors
        try{
            if(tstart[0]) _tstart[0] = PyObject_to_string(tstart[0]);
        } catch (const std::exception& err){ //Catch all errors
            PyErr_SetString(PyExc_RuntimeError,err.what());
            Py_RETURN_NONE;
        }
    }

    //Set the input parameters
    std::vector<std::string> init_params{_mzams[0],_z[0],_spin[0],_snmodel,_tstart[0],
                                         _mzams[1],_z[1],_spin[1],_snmodel,_tstart[1],
                                         _semimajor,_eccentricity,_tend,"all"};

    //Set fake IDs
    size_t ID = 0;
    //Initialise binary and evolve taking into account errors
    Binstar_property *bprop;
    std::string binary_name;
    long binary_rseed = 0;
    double star_tlife[2]={0.,0.};
    try{
        //Initialise the binary
        Binstar b(sevnio->get_sevnio(),init_params,ID,rseed);
        binary_name  = b.get_name();
        binary_rseed = b.get_rseed();
        for (size_t i=0; i<2; i++){
            Star* star = b.getstar(i);
            //Update value
            if (Mcurrent[i] and Mcurrent[i]!=Py_None){
                double _Mass = PyObject_to_double(Mcurrent[i]);
                star->update_from_binary(Mass::ID,_Mass-star->getp(Mass::ID));
            }
            if (MHEcurrent[i] and MHEcurrent[i]!=Py_None){
                double _MHE = PyObject_to_double(MHEcurrent[i]);
                star->update_from_binary(MHE::ID,_MHE-star->getp(MHE::ID));
            }
            if (MCOcurrent[i] and MCOcurrent[i]!=Py_None){
                double _MCO = PyObject_to_double(MCOcurrent[i]);
                star->update_from_binary(MCO::ID,_MCO-star->getp(MCO::ID));
            }
            if (Rcurrent[i] and Rcurrent[i]!=Py_None){
                double _RCO = PyObject_to_double(Rcurrent[i]);
                star->update_from_binary(Radius::ID,_RCO-star->getp(Radius::ID));
            }
        }
        //If just init just initialise to tini
        if (just_init) bprop = new Binstar_property(_init_binary(b));
        else bprop = new Binstar_property(_evolve_binary(b));

        //Estimate tlife post evolution (because maybe meanwhile it changed trakcs)
        star_tlife[0] = b.getstar(0)->get_tphase()[Lookup::Phases::Remnant];
        star_tlife[1] = b.getstar(1)->get_tphase()[Lookup::Phases::Remnant];
    }
    catch (const std::exception& err){ //Catch all errors
        PyErr_SetString(PyExc_RuntimeError,err.what());
        Py_RETURN_NONE;
    }

    //Create the output Dictionary containing the properties
    PyObject* outputProperyDictionary = PyDict_New();
    //Binary parameters
    for (const auto& propID : Binstar_property::bpropertiesID_list){
        auto pname = BinaryProperty::all[propID]->name();
        if (propID==BWorldtime::ID) pname = "Worldtime";
        PyDict_SetItemString(outputProperyDictionary, pname.c_str(), bprop->get(propID));
    }
    //Star parameters
    for (const auto& propID : Star_property::propertiesID_list){
        //Do not store Worldtime
        if (propID!=Worldtime::ID){
            auto pname = Property::all[propID]->name();
            auto pname0 = pname+"_0";
            auto pname1 = pname+"_1";
            PyDict_SetItemString(outputProperyDictionary, pname0.c_str(), bprop->gets(propID,0));
            PyDict_SetItemString(outputProperyDictionary, pname1.c_str(), bprop->gets(propID,1));
        }
    }

    //Create the output Dictionary containing the additional informations
    PyObject* outputAuxiliaryDictionary = PyDict_New();
    PyDict_SetItemString(outputAuxiliaryDictionary,"Log",PyUnicode_FromString(sevnio->get_log().c_str()));
    PyDict_SetItemString(outputAuxiliaryDictionary,"rseed",PyLong_FromLong(binary_rseed));
    PyDict_SetItemString(outputAuxiliaryDictionary,"ID",PyLong_FromSsize_t(ID));
    PyDict_SetItemString(outputAuxiliaryDictionary,"name",PyUnicode_FromString(binary_name.c_str()));
    PyDict_SetItemString(outputAuxiliaryDictionary,"tlife",
                         PyTuple_Pack(2,PyFloat_FromDouble(star_tlife[0]),PyFloat_FromDouble(star_tlife[1])));

    //Reset log
    sevnio->reset_log();
    //Clean bprop
    delete bprop;
    bprop = nullptr;

    PyObject* return_tuple = PyTuple_Pack(2,outputProperyDictionary,outputAuxiliaryDictionary);
    //PyTuple_Pack increments by one the ref count, so just decref it to maintain a value of 1
    Py_DECREF(outputProperyDictionary);
    Py_DECREF(outputAuxiliaryDictionary);

    return return_tuple;

}

///CREATE MODULE
static PyMethodDef callMethods[] = {
        //{"evolve_binary",evolve_binary,METH_NOARGS},
        {"sevnio_initialise",sevnio_initialise,METH_VARARGS},
        {"sevnio_finalise",sevnio_finalise,METH_NOARGS},
        {"sevnio_param",sevnio_param,METH_NOARGS},
        {"evolve_star",reinterpret_cast<PyCFunction>(evolve_star),METH_VARARGS | METH_KEYWORDS},
        {"evolve_binary",reinterpret_cast<PyCFunction>(evolve_binary),METH_VARARGS | METH_KEYWORDS},
        {NULL,NULL,0,NULL}
};

static struct PyModuleDef sevnwrapModule = {
        PyModuleDef_HEAD_INIT,
        "sevnwrap",
        "Test",
        -1,
        callMethods
};

PyMODINIT_FUNC
PyInit_sevnwrap(void) {
    PyObject* m = PyModule_Create(&sevnwrapModule);
    import_array();
    if ( m == NULL ){
        return NULL;
    }


    PyObject* sevnParams;
    PyObject* sevnParamsDescription;
    PyObject* sevnParams_tuple = sevnio_param(NULL,NULL);
    if (!sevnParams_tuple){
        PyErr_Print(); //Print the Error, not need to decref sevnParams_tuple, because GetItem returns a borrowed ref
        Py_DECREF(m);
        return NULL;
    }

    sevnParams    = PyTuple_GetItem(sevnParams_tuple,0);
    //Not needed to increase the sevnParams reference counter, because,
    //when calling AddObject Python already take care of this variable in the module
    if (PyModule_AddObject(m, "sevnParams", sevnParams)<0){
        PyErr_Print(); //Print the Error
        Py_DECREF(sevnParams);
        Py_DECREF(m);
        return NULL;
    }

    sevnParamsDescription    = PyTuple_GetItem(sevnParams_tuple,1);
    if (PyModule_AddObject(m, "sevnParamsDescription", sevnParamsDescription)<0){
        PyErr_Print(); //Print the Error
        Py_DECREF(sevnParamsDescription);
        Py_DECREF(m);
        return NULL;
    }
    //sevnParams_tuple not needed anymore, decref
    Py_XDECREF(sevnParams_tuple);

    ///Add constants
    //G
    PyObject* senvconst_G = PyFloat_FromDouble(utilities::G);
    if (PyModule_AddObject(m, "sevnconst_G", senvconst_G)<0){
        PyErr_Print(); //Print the Error
        Py_DECREF(senvconst_G);
        Py_DECREF(m);
        return NULL;
    }

    //Rsun
    PyObject* sevnconst_Rsun = PyFloat_FromDouble(utilities::Rsun_cgs);
    if (PyModule_AddObject(m, "sevnconst_Rsun", sevnconst_Rsun)<0){
        PyErr_Print(); //Print the Error
        Py_DECREF(sevnconst_Rsun);
        Py_DECREF(m);
        return NULL;
    }



    return m;
}
