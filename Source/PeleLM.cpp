#include <PeleLM.H>

using namespace amrex;

PeleLM::PeleLM() = default;

PeleLM::~PeleLM()
{
   for (int lev = 0; lev <= finest_level; ++lev) {
      ClearLevel(lev);
   }
   prob_parm.reset();
}

PeleLM::LevelData*
PeleLM::getLevelDataPtr(int lev, const PeleLM::TimeStamp &a_time, int useUMac)
{
   AMREX_ASSERT(a_time==AmrOldTime || a_time==AmrNewTime || a_time==AmrHalfTime);
   if ( a_time == AmrOldTime ) { 
      return m_leveldata_old[lev].get();
   } else if ( a_time == AmrNewTime ) {
      return m_leveldata_new[lev].get();
   } else {
      m_leveldata_floating.reset( new LevelData(grids[lev], dmap[lev], *m_factory[lev],
                                  m_incompressible, m_has_divu,
                                  m_nAux, m_nGrowState, m_nGrowMAC));
      Real time = getTime(lev,a_time);
      if (useUMac) {
         // TODO: find a way to get U^{n+1/2} from Umac
         // For now get old time
         Real oldtime = getTime(lev,AmrOldTime);
         fillpatch_velocity(lev, oldtime, m_leveldata_floating->velocity, m_nGrowState);
      } else {
         fillpatch_velocity(lev, time, m_leveldata_floating->velocity, m_nGrowState);
      }
      if (!m_incompressible) {
         fillpatch_density(lev, time, m_leveldata_floating->density, m_nGrowState);
         fillpatch_species(lev, time, m_leveldata_floating->species, m_nGrowState);
         fillpatch_energy(lev, time, m_leveldata_floating->rhoh, m_leveldata_floating->temp, m_nGrowState);
      }
      return m_leveldata_floating.get();
   }
}

PeleLM::LevelDataReact*
PeleLM::getLevelDataReactPtr(int lev)
{
   if (m_do_react) {
      return m_leveldatareact[lev].get();
   } else {
      return nullptr;
   }
}

Vector<MultiFab *>
PeleLM::getVelocityVect(const TimeStamp &a_time) {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->velocity));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->velocity));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getSpeciesVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->species));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->species));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getDensityVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->density));
      }
   } else if ( a_time == AmrNewTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->density));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         Real time = getTime(lev,a_time);
         m_halfTimeDensity[lev].reset( new MultiFab(grids[lev], dmap[lev], 1, m_nGrowState) );
         fillpatch_density(lev,time,*(m_halfTimeDensity[lev].get()),m_nGrowState);
         r.push_back(m_halfTimeDensity[lev].get());
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getTempVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->temp));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->temp));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getDiffusivityVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->diff_cc));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->diff_cc));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getViscosityVect(const TimeStamp &a_time) {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->visc_cc));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->visc_cc));
      }
   }
   return r;
}
