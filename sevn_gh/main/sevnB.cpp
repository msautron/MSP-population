#include <sevn.h>
using sevnstd::SevnLogging;


#include <chrono>
#include <thread>
#include <memory>


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



    //MEMORY DEBUGGING
    //utilities::hardwait("Before loading tables and star",__FILE__,__LINE__);


    //IO class
    IO sevnio(argc, argv);


    //Set the evolve function
    evolve_utility::EvolveFunctor *evolve_function;
    evolve_function = new evolve_utility::EvolveDefault(svlog, true, true);
    //evolve_function = new evolve_utility::EvolveDebug(svlog,true,true);
    //evolve_function = new evolve_utility::EvolveBinaryCompact(svlog,true);
    //evolve_function = new evolve_utility::EvolveDefault(svlog,true,true);
    //evolve_function = new evolve_utility::EvolveBinaryCompact(svlog,true);
    //evolve_function = new evolve_utility::EvolveBLC(svlog,true);
    //evolve_function = new evolve_utility::EvolveBBH(svlog,true);
    //evolve_function = new evolve_utility::EvolveBCX1(svlog,true);
    //evolve_function = new evolve_utility::EvolveW1(svlog,true);
    //evolve_function = new evolve_utility::EvolveBHMS(svlog,true);




    /*
    std::vector<double> v{10,15,16,20};
    utilities::ListGenerator ttest(v);
    std::cout<<ttest.get()<<std::endl;
    ttest.next();
    std::cout<<ttest.get()<<" "<<ttest.empty()<<std::endl;
    ttest.next();
    std::cout<<ttest.get()<<" "<<ttest.empty()<<std::endl;
    ttest.next();
    std::cout<<ttest.get()<<" "<<ttest.empty()<<std::endl;
    ttest.next();
    std::cout<<ttest.empty()<<std::endl;

    utilities::hardwait();
     */

    //MEMORY DEBUGGING
    //utilities::hardwait("After loading tables and star",__FILE__,__LINE__);


    ///Evolve
    //Star clock
    std::chrono::steady_clock::time_point clock_begin = std::chrono::steady_clock::now();
    //Evolve
    int Nchunk = int(sevnio.svpar.get_num("ev_Nchunk"));
    std::vector<Binstar> binaries;
    evolve_utility::chunk_dispatcher(evolve_function,Nchunk,sevnio,binaries,true);
    //Old depracated version
    //evolve_utility::chunk_dispatcher(Nchunk,sevnio,binaries,true,true);

    //End clock
    std::chrono::steady_clock::time_point clock_end = std::chrono::steady_clock::now();


    delete evolve_function;

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
