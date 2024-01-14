/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: German Samolyuk (ORNL) and
                         Mario Pinto (Computational Research Lab, Pune, India)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "compute_heat_flux.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "group.h"
#include "error.h"
#include "kspace.h"
#include "memory.h"

using namespace LAMMPS_NS;

#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

ComputeHeatFlux::ComputeHeatFlux(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{

    bool penflag=false;
    int type;

  //check args, no parital enthalpy
  if (narg < 6) error->all(FLERR,"Illegal compute heat/flux command");

  // check args, yes partial enthalpy
  if(narg>6){
    if(strcmp(arg[6],"penthalpy")==0){
        if ((narg-7)%4==0) penflag=true;  // expect 4 terms per entry
        else error->all(FLERR,"Illegal penthalpy option for compute heat/flux command");
    }
    else error->all(FLERR,"Illegal compute heat/flux command");
  }

   // allocate per-type partial enthalpy terms, with +1 offset
   memory->create(pen, atom->ntypes+1, 3, "heatflux:pen");

   // initialize penthalpy terms to zero
   for (int i=0; i<=atom->ntypes; i++)
     for (int j=0; j<3; j++)
        pen[i][j] = 0.0;

    // overwrite zero values if needed

   if(penflag){
     error->message(FLERR,"Substracting partial enthalpy");

     int ipen=6; // begin args for penthalpy option
     while(ipen<narg-1){
       type=atoi(arg[++ipen]); // atom type

       if(type>atom->ntypes) error->all(FLERR,"Illegal penthalpy option for compute heat/flux command");

       pen[type][0]=atof(arg[++ipen]); // potential
       pen[type][1]=atof(arg[++ipen]); // kinetic
       pen[type][2]=atof(arg[++ipen]); // virial

     }
   }

  vector_flag = 1;
  size_vector = 12;
  extvector = 1;

  // store ke/atom, pe/atom, stress/atom IDs used by heat flux computation
  // insure they are valid for these computations

  int n = strlen(arg[3]) + 1;
  id_ke = new char[n];
  strcpy(id_ke,arg[3]);

  n = strlen(arg[4]) + 1;
  id_pe = new char[n];
  strcpy(id_pe,arg[4]);

  n = strlen(arg[5]) + 1;
  id_stress = new char[n];
  strcpy(id_stress,arg[5]);

  int ike = modify->find_compute(id_ke);
  int ipe = modify->find_compute(id_pe);
  int istress = modify->find_compute(id_stress);


  if (ike < 0 || ipe < 0 || istress < 0)
    error->all(FLERR,"Could not find compute heat/flux compute ID");
  if (strcmp(modify->compute[ike]->style,"ke/atom") != 0)
    error->all(FLERR,"Compute heat/flux compute ID does not compute ke/atom");
  if (modify->compute[ipe]->peatomflag == 0)
    error->all(FLERR,"Compute heat/flux compute ID does not compute pe/atom");
  if (modify->compute[istress]->pressatomflag == 0)
    error->all(FLERR,"Compute heat/flux compute ID does not compute stress/atom");

  vector = new double[12];

}

/* ---------------------------------------------------------------------- */

ComputeHeatFlux::~ComputeHeatFlux()
{
  delete [] id_ke;
  delete [] id_pe;
  delete [] id_stress;
  delete [] vector;
  memory->destroy(pen);
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFlux::init()
{
  // error checks

  int ike = modify->find_compute(id_ke);
  int ipe = modify->find_compute(id_pe);
  int istress = modify->find_compute(id_stress);
  if (ike < 0 || ipe < 0 || istress < 0)
    error->all(FLERR,"Could not find compute heat/flux compute ID");

  c_ke = modify->compute[ike];
  c_pe = modify->compute[ipe];
  c_stress = modify->compute[istress];
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFlux::compute_vector()
{

  invoked_vector = update->ntimestep;

  // invoke 3 computes if they haven't been already

  if (!(c_ke->invoked_flag & INVOKED_PERATOM)) {
    c_ke->compute_peratom();
    c_ke->invoked_flag |= INVOKED_PERATOM;
  }
  if (!(c_pe->invoked_flag & INVOKED_PERATOM)) {
    c_pe->compute_peratom();
    c_pe->invoked_flag |= INVOKED_PERATOM;
  }
  if (!(c_stress->invoked_flag & INVOKED_PERATOM)) {
    c_stress->compute_peratom();
    c_stress->invoked_flag |= INVOKED_PERATOM;
  }

  // heat flux vector = jc[6] + jv[3]
  // jc[6] = convective portion of heat flux = sum_i ke_i*v_i[3], sum pe_i*v_i[3]
  // jv[3] = virial portion of heat flux = sum_i (stress_tensor_i . v_i[3])
  // normalization by volume is not included

  double *ke = c_ke->vector_atom;
  double *pe = c_pe->vector_atom;
  double **stress = c_stress->array_atom;

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int type;

  double jc[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double jv[3] = {0.0,0.0,0.0};

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      type=atom->type[i];

      jc[0] += (pe[i]-pen[type][0])*v[i][0]; //potential
      jc[1] += (pe[i]-pen[type][0])*v[i][1];
      jc[2] += (pe[i]-pen[type][0])*v[i][2];
      jc[3] += (ke[i]-pen[type][1])*v[i][0]; //kinetic
      jc[4] += (ke[i]-pen[type][1])*v[i][1];
      jc[5] += (ke[i]-pen[type][1])*v[i][2];
      jv[0] -= stress[i][0]*v[i][0] + stress[i][3]*v[i][1] + stress[i][4]*v[i][2] - pen[type][2]*v[i][0]; //virial
      jv[1] -= stress[i][3]*v[i][0] + stress[i][1]*v[i][1] + stress[i][5]*v[i][2] - pen[type][2]*v[i][1];
      jv[2] -= stress[i][4]*v[i][0] + stress[i][5]*v[i][1] + stress[i][2]*v[i][2] - pen[type][2]*v[i][2];

    }
  }

// convert jv from stress*volume to energy units via nktv2p factor

  double nktv2p = force->nktv2p;
  jv[0] /= nktv2p;
  jv[1] /= nktv2p;
  jv[2] /= nktv2p;

  // sum across all procs
  // 1st 3 terms are total heat flux
  // 2nd 6 terms are just conductive portion

  double data[12] = {jv[0]+jc[0]+jc[3],jv[1]+jc[1]+jc[4],jv[2]+jc[2]+jc[5],jc[0],jc[1],jc[2],jc[3],jc[4],jc[5],jv[0],jv[1],jv[2]};
  MPI_Allreduce(data,vector,12,MPI_DOUBLE,MPI_SUM,world);
}
