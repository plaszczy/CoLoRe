
#include"Window.hh"

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

#include<string>
#include<vector>
using namespace std;


class Catalog
{
  
public:
  Catalog(string fn,bool rsd=true){
    fitshandle fh;
    fh.open(fn);
    fh.goto_hdu(2);
    //read data arrays
    fh.read_entire_column(2,ra);
    fh.read_entire_column(3,dec);
    fh.read_entire_column(4,z);
    if (rsd) {
      arr<float> dz;
      fh.read_entire_column(5,dz);
      for (size_t i=0;i<dz.size();i++) z[i]+=dz[i];
    }
    cout << "OK: read catalog= " << fn << " with " << z.size() << " entries" << endl;
    fh.close();
      
  }
 

  arr<float> ra,dec,z;
};


class Shell
{

public:

  //constructors
  Shell(const Catalog& c,const Window* w):cat(c),win(w){};

  inline void fillIndex(const size_t& i) {index.push_back(i);}

  vector<size_t> index;
  const Catalog& cat;
  const Window* win;

};


auto main(int argc,char** argv)-> int {

  announce(argv[0]);
  vector<double> zmean={0.15,0.25,0.35,0.45};
  double width=0.05;
  
  bool rsd=true;

  try{
    
    //read catalog
    Catalog cat(argv[1],rsd);
    
    //create shells
    vector<Shell> shell;
    for (size_t i=0;i<zmean.size();i++) 
      shell.push_back(Shell(cat,new UniformWindow(zmean[i]-width,zmean[i]+width)));


    //fill index looping once on the catalog
    for (size_t i=0;i<cat.z.size();i++){
      double z=cat.z[i];
      for (auto &s : shell)
	if (s.win->in(z)) s.fillIndex(i);
    }

    
    for (auto s : shell)
      cout <<" shell [" << s.win->zmin() <<"," <<  s.win->zmax() << "]: " << s.index.size() << " gals" << endl;

    



  }
  catch (PlanckError& e){   
    cerr << e.what()<<endl;
    return 1;
  }
  return 0;

}
