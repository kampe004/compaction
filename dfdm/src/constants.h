#ifndef CONSTANTS_H
#define CONSTANTS_H

/* contains physical constants */

namespace DSM{ 

static const double rho_i=917.; // ice density [kg/m**3]
static const double rho_w=1000.; // water density [kg/m**3]
static const double g=9.81; // earth gravity [m/s**2]
static const double R=8.314; // gas constant [J/mol/K]
static const double T0 = 273.15; // water freezing point [K]
static const double tcair = 0.023; // air thermal conductivity [W/m/K]
static const double tcice = 2.290; // ice thermal conductivity [W/m/K]
static const double cp = 2000.0; // specific heat capacity of ice [J/kg/K] http://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html

static const long sec_in_year=365*24*3600; // amount of seconds in a year
static const long sec_in_day=24*3600; // amount of seconds in a day

static const double fs_dend = 1.0; // fresh snow dendricity 
static const double fs_sphere = 0.5; // fresh snow sphericity
static const double fs_gs   = 0.35; // Brun 1992: when dendricity reaches 0, grain size ranges 0.3-0.4 mm

}

#endif
