#include <PeleLM.H>
#include <pelelm_prob.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   pp.query("P_mean",   prob_parm->P_mean);
   pp.query("inj_start",   prob_parm->inj_start);
   pp.query("inj_dur",   prob_parm->inj_dur);
   pp.query("v_in",  prob_parm->v_in);
   pp.query("D",   prob_parm->D);
   pp.query("Z",   prob_parm->Z);
   pp.query("T_fu",  prob_parm->T_fu);
   pp.query("T_ox",  prob_parm->T_ox);
   pp.query("tau",  prob_parm->tau);
   std::string fu_spec = "";
   std::string fu_ox_spec = "";
   pp.query("fu_spec",   fu_spec);
   pp.query("fu_ox_spec",   fu_ox_spec);
   amrex::Real Y_O2_ox = {0.};
   amrex::Real Y_fu_ox = {0.};
   pp.query("Y_O2_ox",   Y_O2_ox);
   pp.query("Y_fu_ox",   Y_fu_ox);
   
   amrex::Vector<std::string> sname;
   pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(sname);
   amrex::Real Y_pure_fuel[NUM_SPECIES] = {0.0};
   int fu_indx = -1;
   int o2_indx = -1;
   int n2_indx = -1;
   int fu_ox_indx = -1;
   for (int n=0; n<sname.size(); n++)
   {
     if ( sname[n] == fu_spec) fu_indx = n;
     if ( sname[n] == "O2") o2_indx = n;
     if ( sname[n] == "N2") n2_indx = n;
     if ( sname[n] == fu_ox_spec) fu_ox_indx = n;
   }

   if (fu_indx < 0) amrex::Abort("Fuel species not found.");

   Y_pure_fuel[fu_indx] = 1.0;

   prob_parm->Y_ox[o2_indx] = Y_O2_ox;
   prob_parm->Y_ox[fu_ox_indx] = Y_fu_ox;
   prob_parm->Y_ox[n2_indx] = 1.0 - Y_fu_ox - Y_O2_ox;

   for (int n = 0; n < NUM_SPECIES; n++)
   {
       prob_parm->Y_fuel[n] = prob_parm->Z * Y_pure_fuel[n] + (1.-prob_parm->Z) * prob_parm->Y_ox[n];
   }

   CKHBMS(&prob_parm->T_fu, prob_parm->Y_fuel, &prob_parm->H_fuel);
   CKHBMS(&prob_parm->T_ox, prob_parm->Y_ox,   &prob_parm->H_ox);

   auto problo = geom[0].ProbLo();
   auto probhi = geom[0].ProbHi();
   prob_parm->center_xy[0] = 0.5 * (probhi[0] + problo[0]);
   prob_parm->center_xy[1] = 0.5 * (probhi[1] + problo[1]);

#ifdef PELE_USE_TURBINFLOW
   if (pp.countval("turb_file") > 0) {
#if AMREX_SPACEDIM==2
       amrex::Abort("Turbulence inflow unsupported in 2D.");
#endif
       std::string turb_file = "";
       pp.query("turb_file", turb_file);
       amrex::Real turb_scale_loc = 1.0;
       pp.query("turb_scale_loc", turb_scale_loc);
       amrex::Real turb_scale_vel = 1.0;
       pp.query("turb_scale_vel", turb_scale_vel);

       prob_parm->do_turb = true;

       // Hold nose here - required because of dynamically allocated data in tp
       AMREX_ASSERT_WITH_MESSAGE(prob_parm->tp.tph == nullptr,"Can only be one TurbParmHost");
       prob_parm->tp.tph = new TurbParmHost;

       amrex::Vector<amrex::Real> turb_center = {
         {0.5 * (probhi[0] + problo[0]), 0.5 * (probhi[1] + problo[1])}};
       pp.queryarr("turb_center", turb_center);
       AMREX_ASSERT_WITH_MESSAGE(turb_center.size() == 2, "turb_center must have two elements");
       for (int n = 0; n < turb_center.size(); ++n) {
         turb_center[n] *= turb_scale_loc;
       }

       int turb_nplane = 32;
       pp.query("turb_nplane", turb_nplane);
       AMREX_ASSERT(turb_nplane > 0);
       amrex::Real turb_conv_vel = 1;
       pp.query("turb_conv_vel", turb_conv_vel);
       AMREX_ASSERT(turb_conv_vel > 0);

       init_turbinflow(turb_file, turb_scale_loc, turb_scale_vel, turb_center, turb_conv_vel,
                       turb_nplane, prob_parm->tp);
   }
#endif
}
