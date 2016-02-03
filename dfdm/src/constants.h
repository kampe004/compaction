#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Densification { 

static const double rho_i=917.; // ice density (kg/m^3)
static const double rho_w=1000.; // water density
static const double g=9.8; // Gravitational constant (m/s^2)
static const double R=8.314; // Gas Constant (J/mol/K)
static const double T0 = 273.15;

static const double sec_in_year=365.*24*3600; // amount of seconds in a year

static const double tcair = 0.023; // air thermal conductivity [W/m/K]
static const double tcice = 2.290; // ice thermal conductivity [W/m/K]
static const double cp = 2000.0; // specific heat capacity of ice [J/kg/K] http://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html

}

#endif
