/***
 * Created by Giuliano Iorio on 02/09/23.
 * This file contains the class IO_handler that represent a wrapper to the SEVN class IO
 */


#ifndef SEVN_IOWRAP_H
#define SEVN_IOWRAP_H

#include <Python.h>
#include <vector>
#include <string>
#include <stdexcept>


class IO; // Forward declare class IO to avoid multiple definitionand circular import

// Class
class IO_handler{
    /*
     * Singleton class to handle IO from SEVN.
     * This is a Singleton class, so it is build in the way that only one instance of the class can be present
     * at runtime and the constructor is private. The only way to get a pointer to the unique instance is to use
     * initialise() to build the instance and then getInstance() to get the pointer.
     * The instance can be deleted only calling directly the method close, however the two statically dynamically allocated members
     *    can be deleted also with the method reset():
     *      - static IO* senvio_ptr: pointer to the SEVNIO class, also in this case only one instance can be created at time,
     *      therefore it is created and deleted when needed
     *      - static PyObject* _input_dict: Dictionary containing the input parameter to pass to the SEVN IO class,
     *      it is destroyed and re-created accordingly to sevnio_ptr (they are both nullptr or both not nullptr
     */

public:

    //Remove deleting copy constructor
    IO_handler(const IO_handler& obj) = delete;

    ~IO_handler(){
        reset();
    }

    static void initialise(PyObject* input_dict = nullptr){
        if (instace_ptr== nullptr){
            instace_ptr = new IO_handler(input_dict);
        } else {
            PyErr_WarnEx(PyExc_RuntimeWarning,
                         "Calling initialise of an already initialised SEVN instance. We are finalising the object before reinitialising it",
                         1);
            reset();
            //Recall
            instace_ptr = new IO_handler(input_dict);
        }

    }
    static void close(){
        delete instace_ptr;
        instace_ptr = nullptr;
    }

    static IO_handler* getInstance(){
        if (instace_ptr== nullptr){
            throw std::runtime_error("Cannot get an instance of a non initialised SEVN object, please initialise it");
        }
        return instace_ptr;
    }



    inline IO* get_sevnio(){return senvio_ptr;}
    std::string get_log();
    void reset_log();

    static void reset(){
        if (senvio_ptr) delete senvio_ptr;
        senvio_ptr = nullptr;
        //Decref Python Object that is not used
        if (_input_dict) Py_DECREF(_input_dict);
        _input_dict = nullptr;
        //if (instace_ptr) delete instace_ptr;
        //instace_ptr = nullptr;
    }


    std::vector<std::string> get_input_from_dict() const;

    inline std::string IO_param() const;

    static IO_handler* instace_ptr;


private:
    //Make it threadlocal?
    static PyObject* _input_dict;
    static IO* senvio_ptr;

    IO_handler(PyObject* input_dict = nullptr);


};

// Utility functions
IO_handler* check_and_return_IOinstance();

#endif //SEVN_IOWRAP_H
