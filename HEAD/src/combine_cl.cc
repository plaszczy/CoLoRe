
#include "arr.h"
#include"datatypes.h"
#include"string_utils.h"
#include "fitshandle.h"
#include "lsconstants.h"
#include "announce.h"

#include<iostream>
#include<string>
#include<string>
#include<vector>
using namespace std;

using Cl = arr<double>;


auto main(int argc,char** argv)-> int {

PLANCK_DIAGNOSIS_BEGIN

  announce(argv[0]);

 int Ntot= argc-1;
 cout << " Working on" << Ntot << " fits files" <<endl;

 vector<Cl> cls;
 vector<string> colnames; //ensures same counting than icol
 vector< pair<string,string> > keyval;
 arr<int> ell;
 
 bool first=true;
 
   //loop on all fits files
 for (int ifile=1;ifile<argc;ifile++){
   cout << argv[ifile] << endl;
   fitshandle fh;
   fh.open(argv[ifile]);
   fh.goto_hdu(2);
   const int ncols=fh.ncols();
   if (first){
     //colnames
     for (int icol=1;icol<=ncols;icol++) {
       string name=fh.colname(icol);
       colnames.push_back(name);
       cout << icol << " : " << name << endl;
     }
     //keys
     vector<string> keys;
     fh.get_all_keys(keys);
     string val;
     for (auto k : keys) {
       fh.get_key(k,val);
       keyval.push_back(make_pair(k,val));
     }
     //fill ell
     fh.read_entire_column(1,ell);
     first=false;
     //fill first cl
     Cl cl;
     for (int icol=2;icol<=ncols;icol++) {
       fh.read_entire_column(icol,cl);
       cls.push_back(cl);
     }
     for (auto p : keyval){
       cout << p.first << ":" << p.second << endl;
   }
     first=false;
     continue; //next
   }


   // read+coadd
   for (int icol=2;icol<=ncols;icol++) {
     Cl cl;
     fh.read_entire_column(icol,cl);
     for (tsize j=0;j<cls[icol-2].size();j++) cls[icol-2][j]+=cl[j];
   }   
   
 }

 //write summary
 fitshandle fout;
 fout.create("!clmean.fits");
 std::vector<fitscolumn> cols;

 //first is ell
 cols.push_back(fitscolumn(colnames[0], "",1,planckType<double>()));
 for (size_t i=1;i<colnames.size();i++)
   cols.push_back(fitscolumn(colnames[i],"",1,planckType<double>()));
 
 fout.insert_bintab(cols);
 
 int icol=1;
 fout.write_column(icol++,ell);
 for (auto cl : cls)  {
   //mean
   for (tsize j=0;j<cl.size();j++) cl[j]/=Ntot;
   fout.write_column(icol++,cl);
 }
 //add keys from first pass
 for (auto k: keyval) fout.set_key(k.first,k.second);


 //fitshandle fh;
 //fh.open(argv[1]);
 //fout.copy_header(fh);
 

//  vector<string> keys;
//  fh.get_all_keys(keys);

//  for (auto key : keys) {
//    void* value;
//    PDT type=planckType<string>();
//    //cout << "type=" << type << endl;
//    fh.get_key_void(key,value,type);
//    //fout.set_key_void(key,value,type);
//  }


PLANCK_DIAGNOSIS_END
}
