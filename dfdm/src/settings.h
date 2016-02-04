#ifndef SETTINGS_H
#define SETTINGS_H

namespace Densification{

enum class DensificationMethod {Ar10T, Overburden};

typedef struct Settings Settings; // data container to hold user settings
struct Settings {
    /* general */
    bool heat;
    DensificationMethod dm;

    /* overburden parameters */
    double eta0;
    double c5;
    double c6;
};
}

#endif
