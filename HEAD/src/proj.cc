
#include "healpix_map.h"
#include "alm.h"
#include "powspec.h"
#include "healpix_map_fitsio.h"
#include "alm_healpix_tools.h"
#include "arr.h"
#include "fitshandle.h"
#include "powspec_fitsio.h"
#include "announce.h"

#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

auto main(int argc,char** argv)-> int {

  announce(argv[0]);

  
  try{


  }
  catch (PlanckError& e){   
    cerr << e.what()<<endl;
  }
  return 0;

}
