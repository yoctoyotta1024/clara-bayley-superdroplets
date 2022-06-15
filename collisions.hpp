// #ifndef COLLISIONS
// #define COLLISIONS
#include <vector>
#include <random>
#include <algorithm>
#include "superdroplets.hpp"


namespace dlc = dimless_constants;
using namespace std;

#define SDloop(i,nsupers) for(int i=0; i<nsupers; i++)  //for loop over all superdroplets

  
const double prob_jk_const = 1.5e3*(pow(dlc::R0, 3.0))*coll_tstep/iVOL;




int collide_droplets(int nsupers, int nhalf, int scale_p,
               Superdrop* ptr, vector<int> pvec)
{
  Superdrop* p1 = ptr;
  Superdrop* p2 = ptr;
  double phi = 0;
  int gamma = 0;

  double prob, prob_jk;
  double delta_eps, delta_r, delta_m_sol;


  /* Neccesary for random permutation of vector and number generators */
  random_device rd;
  default_random_engine seed(rd());            // for shuffling p vector
  mt19937 gen(rd()); 
  uniform_real_distribution<> dis(0.0, 1.0);   // for generating random phi value 

  /* 0. Randomly shuffle order of indicies to
  droplet pointers in order to generate random pairs */
  shuffle(pvec.begin(), pvec.end(), seed);   


  /*  --- This loop collides each pair of superdroplet ---  */
  for(int i=0; i<nhalf; i++){
   
    /* 1. asign pointers to pair of superdrops
     that will collide for that (p1 -> eps) >= (p2 -> eps) */
    if ((ptr+pvec[2*i]) -> eps > (ptr+pvec[2*i+1]) -> eps){
      p1 = (ptr+pvec[2*i]);
      p2 = (ptr+pvec[2*i+1]);
    } 
    else { // else when (p1 -> eps) = (p2 -> eps)
      p1 = (ptr+pvec[2*i+1]);
      p2 = (ptr+pvec[2*i]); 
    }

    /* 2. Determine probability of coalescence */
    //prob_jk = dis(gen);
    prob_jk = prob_jk_const*(p1->vol() + p2->vol());
    prob = scale_p*max(p1->eps, p2->eps)*prob_jk;     // scaled probability of pair coalescence (p_alpha)

    /* 3. Monte Carlo Step: randomly determine gamma of coalescence */
    phi = dis(gen);                    // random number phi in range [0,1]
    if (phi<(prob-floor(prob))){
      gamma = floor(prob)+1;
    }
    else if(phi>=(prob-floor(prob))){
      gamma = floor(prob);
    }
    int maxgamma = floor((p1 -> eps)/(p2 -> eps));
    gamma = min(gamma, maxgamma);     // gamma factor as per Shima et al. 2009


    /* 3. coalesce particles if gamma != 0 */
    if (gamma != 0){
      //gamma = 1;
      
      /* 3.(option b) if eps1 = gamma*eps2 collide to 
      make twin SDs with same eps, r and mass */
      
      if (p1 -> eps == gamma*(p2 -> eps)){         

        //cout << p1 << " = " << p2 << endl;
        delta_eps = (p2 -> eps)/2; 
        delta_r = pow((p2 -> r), 3.0) + gamma*(pow((p1 -> r), 3.0));
        delta_r = pow(delta_r, 1.0/3);
        delta_m_sol = (p2 -> m_sol) + gamma*(p1 -> m_sol);
        
        p1 -> eps = delta_eps;
        p2 -> eps = (p2 -> eps) - delta_eps;
  
        p1 -> r = delta_r;
        p2 -> r = delta_r; 
  
        p1 -> m_sol = delta_m_sol;
        p2 -> m_sol = delta_m_sol; 

      }

      /* 3.(option a) if eps1 > gamma*eps2 collide to grow 
      drop2 radius and mass via decreasing eps of drop1 SDs*/
      else if( (p1 -> eps) > gamma*(p2 -> eps) ) {       
 
        //cout << p1 << " > with " << p2 << endl;
        delta_eps = (p1 -> eps) - gamma*(p2 -> eps);
        delta_r = pow((p2 -> r), 3.0) + gamma*(pow((p1 -> r), 3.0));
        delta_r = pow(delta_r, 1.0/3);
        delta_m_sol = (p2 -> m_sol) + gamma*(p1 -> m_sol);

        p1 -> eps = delta_eps;

        p2 -> r = delta_r;
        p2 -> m_sol = delta_m_sol;

      }

        /* 3.(option "c") this should not be an option, throw error warning */
      else{
        cout << "COLLLISION WARNING!! something undefined occured!!" << endl;
        cout << (p1 -> eps) << " < " << gamma*(p2 -> eps)  << " ?" << endl;
      }
     
        
    } // end coalescence


  } // end collisions of all SD pairs
  




  return 0;
}














// #endif // COLLISIONS