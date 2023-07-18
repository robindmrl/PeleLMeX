#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void pelelm_initdata(int i, int j, int k,
                     int is_incompressible,
                     amrex::Array4<amrex::Real> const& state,
                     amrex::Array4<amrex::Real> const& aux,
                     amrex::GeometryData const& geomdata,
                     ProbParm const& prob_parm,
                     pele::physics::PMF::PmfData::DataContainer const * pmf_data)
                     
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

    AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
                 const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    auto eos = pele::physics::PhysicsType::eos();
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real u[3] = {0.0};
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    //amrex::Real y1, y2;
    
    amrex::Real P_cgs = prob_parm.P_mean * 10.0;
    
    // Limit for the 2 zones of the coflow => [ fuel | oxydizer | air ]
    amrex::Real split_fo = prob_lo[0] + prob_parm.fuel_ox_split;
    amrex::Real split_oa = prob_lo[0] + prob_parm.ox_air_split;

    // We set temp, velocity and mass fractions for the 3 zones; then density, enthalpy and mass.
    if (x < split_fo and y <= prob_parm.ymax) // fuel
    {
	/*
	Chercher pour popuvoir changer le fuel et le dilaunt
	*/
	for (int n = 0; n < NUM_SPECIES; n++){
	    massfrac[n] = 0;}
		    
	// mass fractions
	massfrac[CH4_ID] = prob_parm.dilution;
	massfrac[N2_ID] = 1.0-prob_parm.dilution;	
	    
        //Temp
        state(i,j,k,TEMP) = prob_parm.T_fu;
	    
	//Velocity Inflow
	AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
		     state(i,j,k,VELY) = prob_parm.V_fu;,
		     state(i,j,k,VELZ) = 0.0); 
    }
    else if (x < split_oa and y <= prob_parm.ymax) // oxydizer
    {
    	for (int n = 0; n < NUM_SPECIES; n++){
	    massfrac[n] = 0;}
		    
	// mass fractions
	massfrac[O2_ID] = 0.233;
	massfrac[N2_ID] = 0.767;	
	    
        //Temp
        state(i,j,k,TEMP) = prob_parm.T_ox;
	    
	//Velocity Inflow
	AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
		     state(i,j,k,VELY) = prob_parm.V_ox;,
		     state(i,j,k,VELZ) = 0.0); 
    }
    else // air 
    {
    	for (int n = 0; n < NUM_SPECIES; n++){
	    massfrac[n] = 0;}
		    
	// mass fractions
	massfrac[O2_ID] = 0.233;
	massfrac[N2_ID] = 0.767;	
	    
        //Temp
        state(i,j,k,TEMP) = prob_parm.T_air;
	    
	//Velocity Inflow
	AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
		     state(i,j,k,VELY) = 0.0;,
		     state(i,j,k,VELZ) = 0.0);
    }
    
    //The ignition
    if (prob_parm.do_ignition==1)
    {
        amrex::Real xc = prob_parm.fuel_ox_split;
        amrex::Real yc = prob_parm.ymax;
        amrex::Real rad = std::sqrt((y-yc)*(y-yc)+(x-xc)*(x-xc));
        
        amrex::Real Zst = 0.0; //stoechiometry for the given dilution of CH4 into N2
        amrex::Real X=0.0; //for calculus purpose
        amrex::Real dilution_vol = 0.0; //volumic dilution (the input is mass dilution)
        amrex::Real mtot = 0.0;
        
        if (rad <= prob_parm.ign_radius)
        {
            state(i,j,k,TEMP) = prob_parm.ign_T;
            
            X=(16*(1-prob_parm.dilution))/(28*prob_parm.dilution);
            dilution_vol = 1/(X+1);
            Zst = (4.0/0.42)*dilution_vol;
            mtot = 16*prob_parm.ign_rich*dilution_vol + 28*(0.79*Zst+prob_parm.ign_rich*(1-dilution_vol)) + 32*0.21*Zst;
            
            massfrac[O2_ID] = 32*0.21*Zst/mtot; //21.9% pour richesse = 1
            massfrac[N2_ID] = 28*(0.79*Zst+prob_parm.ign_rich*(1-dilution_vol))/mtot; //72.6% pour richesse = 1
            massfrac[CH4_ID] = 16*prob_parm.ign_rich*dilution_vol/mtot; //5.5% pour richesse = 1
        }
    }
    
     // Density
     amrex::Real rho_cgs = 0.0;
     eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
     state(i,j,k,DENSITY) = rho_cgs * 1.0e3;

     // Enthalpy
     amrex::Real h_cgs = 0.0;
     eos.TY2H(state(i,j,k,TEMP), massfrac, h_cgs);
     state(i,j,k,RHOH) = h_cgs * 1.0e-4 * state(i,j,k,DENSITY);

     // Species mass
     for (int n = 0; n < NUM_SPECIES; n++) {
        state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);}
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const int m_nAux,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
   const amrex::Real* prob_lo = geomdata.ProbLo();
   const amrex::Real* prob_hi = geomdata.ProbHi();

   AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
                const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

   amrex::Real split_fo = prob_lo[0] + prob_parm.fuel_ox_split;
   amrex::Real split_oa = prob_lo[0] + prob_parm.ox_air_split;
    
   amrex::Real mass_frac_air[NUM_SPECIES] = {0.0};
   amrex::Real mass_frac_ox[NUM_SPECIES] = {0.0};
   amrex::Real mass_frac_fuel[NUM_SPECIES] = {0.0};
   
   amrex::Real rho_cgs, P_cgs, RhoH_temp;

   auto eos = pele::physics::PhysicsType::eos();

   for (int n = 0; n < NVAR; n++){
      s_ext[n] = 0.0;
   }

// Get zone
   int zone = -1;
   if ( std::abs(x[0]-prob_lo[0])< prob_parm.fuel_ox_split ) 
   {
      zone = 0; // fuel
      mass_frac_fuel[CH4_ID] = prob_parm.dilution;
      mass_frac_fuel[N2_ID] = 1.0-prob_parm.dilution;
   } 
   else if ( std::abs(x[0]-prob_lo[0])< prob_parm.ox_air_split )
   {
      zone = 1; // oxydizer
      mass_frac_ox[O2_ID]=0.233;
      mass_frac_ox[N2_ID]=0.767;
   }
   else 
   {
      zone = 2; // air
      mass_frac_air[O2_ID]=0.233;
      mass_frac_air[N2_ID]=0.767;
   }

   P_cgs = prob_parm.P_mean * 10.0;

   //Inflow velocity
   /*
   Essayer de mettre un profil hyperbolique (voir Hugo)
   */
   if (idir == 1 && sgn == 1) 
   {
      s_ext[VELX] = 0.0;
      
      if ( zone == 0 ) 
      {
         s_ext[VELY] = prob_parm.V_fu;
         s_ext[TEMP] = prob_parm.T_fu;
         
         eos.PYT2R(P_cgs, mass_frac_fuel, s_ext[TEMP], rho_cgs);
         s_ext[DENSITY] = rho_cgs * 1.0e3;

         eos.TY2H(s_ext[TEMP], mass_frac_fuel, RhoH_temp);
         s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

         for (int n = 0; n < NUM_SPECIES; n++) {
            s_ext[FIRSTSPEC+n] = mass_frac_fuel[n] * s_ext[DENSITY];}
      } 
      else if ( zone == 1 ) 
      {
         s_ext[VELY] = prob_parm.V_ox;
         s_ext[TEMP] = prob_parm.T_ox;
         
         eos.PYT2R(P_cgs, mass_frac_ox, s_ext[TEMP], rho_cgs);
         s_ext[DENSITY] = rho_cgs * 1.0e3;

         eos.TY2H(s_ext[TEMP], mass_frac_ox, RhoH_temp);
         s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

         for (int n = 0; n < NUM_SPECIES; n++) {
            s_ext[FIRSTSPEC+n] = mass_frac_ox[n] * s_ext[DENSITY];}
      }
      else 
      {
         s_ext[VELY] = 0.0;
         s_ext[TEMP] = prob_parm.T_air;
         
         eos.PYT2R(P_cgs, mass_frac_air, s_ext[TEMP], rho_cgs);
         s_ext[DENSITY] = rho_cgs * 1.0e3;

         eos.TY2H(s_ext[TEMP], mass_frac_air, RhoH_temp);
         s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

         for (int n = 0; n < NUM_SPECIES; n++) {
            s_ext[FIRSTSPEC+n] = mass_frac_air[n] * s_ext[DENSITY];}
      }      
     }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& domainBox,
           const int  dir,
           const int  beta_comp,
           const int  nComp)
{
     amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);
    // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
    // otherwise this routine could be called for other face diffusivity (Temp, velocity, ...)
}
#endif
