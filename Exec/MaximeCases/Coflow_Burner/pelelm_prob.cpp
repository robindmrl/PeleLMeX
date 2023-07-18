#include <PeleLM.H>
#include <AMReX_ParmParse.H>
#include <PMFData.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   std::string type;
   pp.query("P_mean", PeleLM::prob_parm->P_mean);
   pp.query("standoff", PeleLM::prob_parm->standoff);
   pp.query("pertmag",  PeleLM::prob_parm->pertmag);

   
   pp.query("fuel_ox_split", PeleLM::prob_parm->fuel_ox_split);
   pp.query("ox_air_split", PeleLM::prob_parm->ox_air_split);
   pp.query("ymax",PeleLM::prob_parm->ymax);
   pp.query("V_fu",PeleLM::prob_parm->V_fu);
   pp.query("V_ox",PeleLM::prob_parm->V_ox);
   pp.query("T_fu",PeleLM::prob_parm->T_fu);
   pp.query("T_ox",PeleLM::prob_parm->T_ox);
   pp.query("T_air",PeleLM::prob_parm->T_air);
   pp.query("dilution",PeleLM::prob_parm->dilution);
   pp.query("do_ignition",PeleLM::prob_parm->do_ignition);
   pp.query("ign_radius",PeleLM::prob_parm->ign_radius);
   pp.query("ign_T",PeleLM::prob_parm->ign_T);
   pp.query("ign_rich",PeleLM::prob_parm->ign_rich);

}
