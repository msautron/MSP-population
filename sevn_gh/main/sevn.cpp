#include <sevn.h>
using sevnstd::SevnLogging;

#include <chrono>
#include <thread>
using sevnstd::SevnLogging;

int main (int argc, char **argv){

    //TODO should print out a help function with the option -h
    //A  temporary workaround
    if (argc==1){
        IO sevnio;
        //Print used params
        std::cout<<sevnio.svpar;
        //Print SEVN version info
        std::cout<<SEVNinfo::get_full_info()<<std::endl;
        return EXIT_SUCCESS;
    }

    SevnLogging svlog; //TO initialise the debug level

    //IO class
    IO sevnio(argc, argv);
    //Set the evolve function
    thread_local evolve_utility::EvolveFunctor* evolve_function;
    evolve_function = new evolve_utility::EvolveDefault(svlog,true);



    ///Evolve
    //Star clock
    std::chrono::steady_clock::time_point clock_begin = std::chrono::steady_clock::now();
    //Evolve
    int Nchunk = int(sevnio.svpar.get_num("ev_Nchunk"));
    std::vector<Star> stars;
    evolve_utility::chunk_dispatcher(evolve_function,Nchunk,sevnio,stars,true);
    //Old depreacated version
    //evolve_utility::chunk_dispatcher(Nchunk,sevnio,stars,true,true);
    //End clock
    std::chrono::steady_clock::time_point clock_end = std::chrono::steady_clock::now();


    ///Final
    //Flush
    std::cerr<<std::flush;
    std::cout<<std::flush;
    //Print info
    std::cout<<"Evolution Done"<<std::endl;
    long Total_time = std::chrono::duration_cast<std::chrono::milliseconds>(clock_end - clock_begin).count();
    std::cout<<"Total time: "<<Total_time/1e3<<" s"<<std::endl;
    std::cout<<"Total time: "<<Total_time/(60*1e3)<<" m"<<std::endl;
    std::cout<<"Total time: "<<Total_time/(3600*1e3)<<" h"<<std::endl;
    //Print used params
    std::cout<<sevnio.svpar;
    sevnio.print_params();

    //Print SEVN version info
    std::cout<<SEVNinfo::get_full_info()<<std::endl;

    return 0;

}

