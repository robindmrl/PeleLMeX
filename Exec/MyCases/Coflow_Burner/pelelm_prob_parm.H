#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
   amrex::Real P_mean   = 101325.0_rt;
   amrex::Real standoff = 0.0_rt;
   amrex::Real pertmag  = 0.0004_rt;
   
   amrex::Real fuel_ox_split = 0.0; // distance from axis of fuel/ox limit = radius of fuel flow
   amrex::Real ox_air_split = 0.0; // distance from axis of ox/fuel limit = radius of ox flow
   
   amrex::Real ymax = 0.0; // ymax of inflow in the initial condition
   
   amrex::Real V_fu = 0.0; //fuel inflow velocity
   amrex::Real V_ox = 0.0; // oxydizer inflow velocity
   // air inflow velocity is always set up to 0.0 in our case
   
   amrex::Real T_fu = 298.0; // fuel inflow temp
   amrex::Real T_ox = 298.0; // oxidyzer inflow temp
   amrex::Real T_air = 298.0; // air inflow temp
   
   amrex::Real dilution = 1.0; //part of Ch4 in the fuel tank (CH4+N2)
   
   int do_ignition = 1;
   amrex::Real ign_radius = 0.0; //Radius of the hot bubble for the ignition
   amrex::Real ign_T = 1000.0; //Temp of the hot bubble at the start
   amrex::Real ign_rich = 0.0; //richness of mixture in the bubble : ign_rich*fuel + Zst*Oxydizer = ...
};
#endif
