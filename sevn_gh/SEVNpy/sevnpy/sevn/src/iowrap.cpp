/***
 * Created by Giuliano Iorio on 02/09/23.
 * This file contains the class IO_handler that represent a wrapper to the SEVN class IO
 */

#include "iowrap.h"
#include <Python.h>
#include <IO.h>
#include <vector>
#include <string>
#include "wraputils.h"



IO* IO_handler::senvio_ptr = nullptr;
IO_handler* IO_handler::instace_ptr = nullptr;
PyObject* IO_handler::_input_dict = nullptr;
IO_handler::IO_handler(PyObject* input_dict){

    /** CHeck the input **/
    //is not a nullptr and it is a Dictionary: OK!
    if(input_dict and PyDict_Check(input_dict)){
        _input_dict = input_dict;
        //We will mantain input_dict in the class, so INCREF THE COUNTER
        Py_INCREF(_input_dict);
    }
        //it is None or a nullptr: OK! just mantein _input_dict in the nullptr status
    else if (input_dict and (input_dict==Py_None or input_dict== nullptr)){
        _input_dict = nullptr;
    }
        //is not a nullptr and it is not a Dictionary: NO, throw an error!
    else if(input_dict){
        throw std::runtime_error("Input in IO_handler is not a Dictionary");
    }

    /** Initialise IO class **/
    if (senvio_ptr== nullptr){
        //Initialise
        senvio_ptr = new IO();

        //Get the string input
        auto ioinput = get_input_from_dict();

        //Now define the char array used in the IO load taking the input from the ioinput vector
        //@TODO get_input_from_dict could return directy a charArray
        int nstring = ioinput.size();
        auto charArray = new char*[nstring];
        for (int i =0; i<  nstring;i++){
            charArray[i] = new char[ioinput[i].size() +1];
            strcpy(charArray[i], ioinput[i].c_str());
        }

        //Try to load or in case clean everything and throw
        try{
            senvio_ptr->load(nstring,charArray);
            //Delete the dynamically allocated array of char
            for (int i=0; i < nstring; i++){
                delete[] charArray[i];
            }
            delete[] charArray;
        }
        catch (const std::exception& err) {
            //Delete the dynamically allocated array of char
            for (int i=0; i < nstring; i++){
                delete[] charArray[i];
            }
            delete[] charArray;
            //Finalise
            reset();
            throw;
        }


    }



}
std::string IO_handler::get_log(){return senvio_ptr->get_logstring();}
void IO_handler::reset_log(){senvio_ptr->reset_log();}
std::string IO_handler::IO_param() const{
    std::string output;
    if (senvio_ptr) output = senvio_ptr->svpar.print();
    return output;
}

std::vector<std::string> IO_handler::get_input_from_dict() const {

    std::vector<std::string> output = {"SEVN",};
    //If _input_dict is nullptr just return the 1 element vector
    if (_input_dict==nullptr) return output;

    //If _input_dict is not null takes the values
    //Key, value objects, they are borrowed, not need to INCREF-DECREF (https://docs.python.org/3/c-api/dict.html)
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    //Iterate over pos,key pairs
    while (PyDict_Next(_input_dict, &pos, &key, &value)){
        auto key_str   = PyObject_to_string(key);
        auto value_str = PyObject_to_string(value);
        output.push_back(std::string{"-"}+key_str);
        output.push_back(value_str);
    }


    return output;
}

IO_handler* check_and_return_IOinstance(){

    IO_handler* sevnio;
    try{
        sevnio = IO_handler::getInstance();
        if (sevnio== nullptr) throw std::runtime_error("SEVN has not been initiliased, please run sevnio_initialise().");
    }
    catch (const std::exception& err){ //Catch all errors
        PyErr_SetString(PyExc_RuntimeError,err.what());
        return nullptr;
    }


    return sevnio;
}