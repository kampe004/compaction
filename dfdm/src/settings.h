#ifndef SETTINGS_H
#define SETTINGS_H

namespace Densification{

enum class DensificationMethod {Ar10T, Overburden};

typedef struct Settings Settings; // data container to hold user settings
struct Settings {

    /* physics */
    bool heat;
    DensificationMethod dm;

    /* overburden parameters */
    double eta0;
    double c5;
    double c6;

    /* meteorological forcing */
    int forcing_dt; 
    std::string f_acc; 
    std::string f_wind10m;
    std::string f_tskin;
};
}

#endif
