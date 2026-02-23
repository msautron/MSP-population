//
// Created by iorio on 5/21/24.
//

#ifndef SEVN_WINDACCRETION_H
#define SEVN_WINDACCRETION_H

#include <IO.h>
#include <Processes.h>
#include <binstar.h>
#include <star.h>

//Forward declaration
class Star;
class Binstar;


class Windaccretion : public MaccretionProcess {

public:
    Windaccretion(_UNUSED IO* _io = nullptr, bool reg = true) {
        if (reg) {
            //First register to include this Property in static array
            Register(this, &ID, name());
        }
        betaw = 0.125;
        alphaw = 1.5;
        muw = 1.0;

        if (_io!= nullptr){
            alphaw = _io->svpar.get_num("w_alpha");
            betaw  = _io->svpar.get_num("w_beta");
        }
    }

    static size_t ID;
    static Windaccretion _windaccretion;
    inline std::string name() override { return "Windaccretion"; }

    Windaccretion* Instance(_UNUSED IO* _io) override;
    virtual Windaccretion* instance() { return nullptr; }


    ///----------

    static constexpr int NO_WINDACCRETION = 0;
    virtual double DA(_UNUSED Binstar* b, _UNUSED int procID) { return 0; }
    virtual double DE(_UNUSED Binstar* b, _UNUSED int procID) { return 0; }
    inline bool is_process_ongoing() const override { return false; }

    ///-------------

protected:

    static std::map<std::string, Windaccretion*>& GetStaticMap() {
        static std::map<std::string, Windaccretion*> _locmap;
        return _locmap;
    }
    static std::vector<int>& GetUsed() {
        static std::vector<int> _used;
        return _used;
    }
    static void Register_specific(Windaccretion* ptr) {
        GetStaticMap().insert(std::make_pair(ptr->name(), ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size() - 1] = 0;
    }

    double betaw; /*!< Wind escape velocity parameter (Eq. 9, Hurley+02) */
    double alphaw; /*!< Bondi-Hoyle accretion parameter (Eq. 6, Hurley+02) */
    double muw; /*!< Angular momentum  transfer efficiency  (Eq. 11, Hurley+02) */

};

class WindaccretionHurley : public Windaccretion {

public:

    WindaccretionHurley(_UNUSED IO* _io = nullptr, bool reg = true) : Windaccretion(nullptr, false) {
        if (reg) {
            Register_specific(this);
        }
    }

    inline std::string name() override { return "hurley"; }
    static WindaccretionHurley _windaccretionhurley;

    WindaccretionHurley* instance() override {
        return (new WindaccretionHurley(nullptr, false));
    }

    int evolve(Binstar* binstar) override;

    double estimate_accreted_mass(_UNUSED  double DM,  _UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *binstar) const  override;
    double DA(_UNUSED Binstar* b, _UNUSED int procID) override;
    double DE(_UNUSED Binstar* b, _UNUSED int procID) override;
    inline bool is_process_ongoing() const override { return true; }

protected:

    int accrete_mass(Binstar* binstar) override;


};

class WindaccretionDisabled : public Windaccretion {

public:

    WindaccretionDisabled(_UNUSED IO* _io = nullptr, bool reg = true) : Windaccretion(nullptr, false) {
        if (reg) {
            //Second register to handle the different options
            Register_specific(this);
        }
    }

    inline std::string name() override { return "disabled"; }
    static WindaccretionDisabled _windaccretiondisabled;
    WindaccretionDisabled* instance() override {
        return (new WindaccretionDisabled(nullptr, false));
    }
    int evolve(_UNUSED Binstar* binstar) override;

    inline bool is_process_ongoing() const override { return false; }
};



#endif //SEVN_WINDACCRETION_H
