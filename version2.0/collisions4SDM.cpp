// Author: Clara Bayley
// File: collisions4SDM.cpp
/* Functionality for modelling
  collisions in superdroplet model */


#include "collisions4SDM.hpp"




static Superdrop& get_dropletref(const string whichdrop,
                Superdrop &dropA, Superdrop &dropB)
  /* compare dropA.eps with dropB.eps and return 
  either drop1 or drop2 such that drop1.eps is always >
  drop2.eps */
{    
  
  if (dropA.eps > dropB.eps)
  {
    if (whichdrop == "drop1") { return dropA; }
    else { return dropB; }
  }

  else  
  {
    if (whichdrop == "drop1") { return dropB; }
    else { return dropA; }
  }

}





void collide_droplets(const int nsupers, const int nhalf, vector<int> pvec,
       Superdrop (&superdrops_arr)[init::NSUPERS])
{
  double prob, prob_jk;
  double delta_eps, delta_r, delta_m_sol;
  double phi = 0;
  int gamma = 0;
  
  static const double prob_jk_const = 1.5e3*(pow(dlc::R0, 3.0))*init::COLL_TSTEP/init::DROPVOL;
  const int scale_p = nsupers*(nsupers-1)/(2*nhalf);


  /* Neccesary for random permutation of vector and number generators */
  random_device rd;
  default_random_engine seed(rd());            // for shuffling p vector
  mt19937 gen(rd()); 
  uniform_real_distribution<> dis(0.0, 1.0);   // for generating random phi value 

  /* Randomly shuffle order of indicies to
  droplet pointers in order to generate random pairs */
  shuffle(pvec.begin(), pvec.end(), seed);   


  /*  --- This loop collides each pair of superdroplets ---  */
  for(int i=0; i<nhalf; i++)
  {
   
    /* 1. asign pointers to pair of superdrops
     that will collide for that (drop1.eps) >= (drop2.eps) */
    Superdrop& drop1 = get_dropletref("drop1", superdrops_arr[pvec[2*i]],
                                  superdrops_arr[pvec[2*i+1]]);
    Superdrop& drop2 = get_dropletref("drop2", superdrops_arr[pvec[2*i]],
                                  superdrops_arr[pvec[2*i+1]]);

    /* 2. Determine probability of coalescence */
    //prob_jk = dis(gen);
    prob_jk = prob_jk_const*(drop1.vol() + drop2.vol());
    prob = scale_p*max(drop1.eps, drop2.eps)*prob_jk;     // scaled probability of pair coalescence (p_alpha)

    /* 3. Monte Carlo Step: randomly determine gamma of coalescence */
    phi = dis(gen);                    // random number phi in range [0,1]
    if (phi<(prob-floor(prob))){
      gamma = floor(prob)+1;
    }
    else if(phi>=(prob-floor(prob))){
      gamma = floor(prob);
    }
    int maxgamma = floor((drop1.eps)/(drop2.eps));
    gamma = min(gamma, maxgamma);     // gamma factor as per Shima et al. 2009


    /* 3. coalesce particles if gamma != 0 */
    if (gamma != 0){
      //gamma = 1;
      
      /* 3.(option b) if eps1 = gamma*eps2 collide to 
      make twin SDs with same eps, r and mass */
      
      if (drop1.eps == gamma*(drop2.eps)){         

        //cout << p1 << " = " << p2 << endl;
        delta_eps = (drop2.eps)/2; 
        delta_r = pow((drop2. r), 3.0) + gamma*(pow((drop1.r), 3.0));
        delta_r = pow(delta_r, 1.0/3);
        delta_m_sol = (drop2. m_sol) + gamma*(drop1.m_sol);
        
        drop1.eps = delta_eps;
        drop2. eps = (drop2. eps) - delta_eps;
  
        drop1.r = delta_r;
        drop2.r = delta_r; 
        
        drop1.m_sol = delta_m_sol;
        drop2.m_sol = delta_m_sol; 

      }

      /* 3.(option a) if eps1 > gamma*eps2 collide to grow 
      drop2 radius and mass via decreasing eps of drop1 SDs*/
      else if( (drop1.eps) > gamma*(drop2.eps) ) {       
 
        //cout << p1 << " > with " << p2 << endl;
        delta_eps = (drop1.eps) - gamma*(drop2.eps);
        delta_r = pow((drop2.r), 3.0) + gamma*(pow((drop1.r), 3.0));
        delta_r = pow(delta_r, 1.0/3);
        delta_m_sol = (drop2.m_sol) + gamma*(drop1.m_sol);

        drop1.eps = delta_eps;

        drop2.r = delta_r;
        drop2.m_sol = delta_m_sol;

      }

        /* 3.(option "c") this should not be an option, throw error warning */
      else{
        cout << "COLLLISION WARNING!! something undefined occured!!" << endl;
        cout << (drop1.eps) << " < " << gamma*(drop2.eps)  << " ?" << endl;
      }
     
        
    } // end coalescence


  } // end collisions of all SD pairs
  

}

