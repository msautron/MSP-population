//
// Created by iorio on 1/27/23.
//

#ifndef SEVN_SEVN_VERSION_H
#define SEVN_SEVN_VERSION_H

#include <string>

struct SEVNinfo{

    static constexpr unsigned int VERSION_MAJOR     =  2;
    static constexpr unsigned int VERSION_MINOR     =  10;
    static constexpr unsigned int VERSION_PATCH     =  0;
    static constexpr auto VERSION    = "2.10.0";
    static constexpr auto DEVLINE    = "zelaous_redgiant";

    //GIT BRANCH Name
    static constexpr auto GIT_BRANCH                = "NOT AVAILABLE";

    //GIT Local HEAD
    //These are the info about the local repo not yet pushed to the remote branch
    static constexpr auto GIT_SHA                                = "NOT AVAILABLE";
    static constexpr auto GIT_SHATIME                            = "NOT AVAILABLE";
    static constexpr unsigned long int GIT_SHATIMESTAMP          = 0;
	static constexpr unsigned int GIT_COUNTER                    = 0;

	//GIT Local remote HEAD
	//These are the info about the last updates on the  remote branch
    static constexpr auto GIT_SHA_REMOTE                                = "NOT AVAILABLE";
    static constexpr auto GIT_SHATIME_REMOTE                            = "NOT AVAILABLE";
    static constexpr unsigned long int GIT_SHATIMESTAMP_REMOTE          = 0;
	static constexpr unsigned int GIT_COUNTER_REMOTE                    = 0;

    static std::string get_full_info(){

        std::string full_info;
        full_info+="*********************\n";
        full_info+="SEVN version info:\n";
        full_info+="---------------------\n";
        full_info+=" SEVN DEVELOPMENT LINE:" + std::string(SEVNinfo::DEVLINE);
        full_info+="\n SEVN VERSION:" + std::string(SEVNinfo::VERSION);
        full_info+="\n SEVN BRANCH:" + std::string(SEVNinfo::GIT_BRANCH);
        full_info+="\n---------------------\n";
        full_info+="Git local info:";
        full_info+="\n COMMITS COUNTER:" + std::to_string(SEVNinfo::GIT_COUNTER);
        full_info+="\n LAST COMMIT HASH:" + std::string(SEVNinfo::GIT_SHA);
        full_info+="\n LAST COMMIT DATE:" + std::string(SEVNinfo::GIT_SHATIME)  + " (TIMESTAMP: "+std::to_string(SEVNinfo::GIT_SHATIMESTAMP)+")";
        full_info+="\nGit remote info:";
        full_info+="\n COMMITS COUNTER:" + std::to_string(SEVNinfo::GIT_COUNTER_REMOTE);
        full_info+="\n LAST COMMIT HASH:" + std::string(SEVNinfo::GIT_SHA_REMOTE);
        full_info+="\n LAST COMMIT DATE:" + std::string(SEVNinfo::GIT_SHATIME_REMOTE)  + " (TIMESTAMP: "+std::to_string(SEVNinfo::GIT_SHATIMESTAMP_REMOTE)+")";
        full_info+="\n*********************\n";

        return full_info;
    }
};

#endif //SEVN_SEVN_VERSION_H
