//--------------------------------------------------------------------------
//
// Description:
// 	class Timer :
//
//
// History (add to end):
//	creation:   ven mar 10 19:43:01 CET 2006 
//
//------------------------------------------------------------------------

#ifndef Timer_hh
#define Timer_hh

#include <time.h>      // Needed for unix implementation of timer
#include <sys/time.h>  // Needed for unix implementation of timer
#include <sys/types.h> // Needed for unix implementation of timer
#include <sys/times.h> // Needed for unix implementation of timer

#include<iostream>
#include<string>

class Timer
{

public:

  //constructors
  Timer():t0(clock()),T0(getUniversalTime()){tlast=t0; Tlast=T0;};

  inline float totalCPU(){
    clock_t t=clock();
    return float(t-t0)/CLOCKS_PER_SEC;
  }

  inline float partialCPU(){
    clock_t t=clock();
    float dt=float(t-tlast)/CLOCKS_PER_SEC;
    tlast=t;
    return dt;
  }

  inline float total(){
    double T=getUniversalTime();
    return T-T0;
  }

  inline float partial(){
    double T=getUniversalTime();
    double dt=T-Tlast;
    Tlast=T;
    return dt;
  }


  // destructor
  ~Timer(){};

  static double getUniversalTime() {
    timeval CurrentTime;
    static double scale = 1.0/1000000.0;
    gettimeofday( &CurrentTime, NULL );
    return CurrentTime.tv_sec + scale * CurrentTime.tv_usec;
  }


private:
  clock_t t0;
  clock_t tlast;
  double T0,Tlast;
  
};

inline std::ostream& operator<<(std::ostream& o,Timer& t){
  o << " [time: partial="<< t.partial() << " s /total="<< t.total() << " s]";
  return o;
}

#endif

