#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

#include <cmath>

#include <vector>
#include <map>

#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <sstream>
#include <iterator>

#include "common.hpp"

bool file_exists(const std::string& filename)
{
  struct stat buf;
  if (stat(filename.c_str(), &buf) != -1)
    return true;
                    
  return false;
}

int fill_matrix (Matrix & x, Matrix & y, std::vector<float> & logp_v, 
    std::vector<std::vector<int> > & finger_v)
{

  int nobj = (int) logp_v.size();
  int descdim = (int) finger_v.begin()->size();

  lsAllocMatrix (&x, nobj, descdim);
  lsAllocMatrix (&y, nobj, 1);
  for (int i=0; i<nobj; ++i)
  {
    y.data[i][0] = logp_v[i];
    for (int j=0; j<(int)finger_v[i].size(); ++j)
      x.data[i][j] = finger_v[i][j];
  }
 
  return nobj;
}

int build_pls (Matrix & x, Matrix & y, int num_of_components, 
    float * sdec, double * r2, LV * plscoeff)
{
  float vx, vxac, sy;
  int comp, i;
  
  TLV tlv;
  
  vxac  = 0.00;
  
  lsCheckMissing (&x, -99.00);
  lsCheckMissing (&y, -99.00);
  
  lsAutoscaleMatrix (&y);
  lsAutoscaleMatrix (&x);
  
  lsAllocTLV (&tlv, &x, &y);
  lsInitTLV (&tlv, &x, &y);
  
  for (comp=0; comp<=num_of_components; comp++) 
  {
    lsAllocLV (&plscoeff[comp], &x, &y);
    
    lsExtractLV  (&x, &y, &plscoeff[comp]);
    lsDeflatePLS (&x, &y, &plscoeff[comp], &tlv, (comp+1));
  
    vx = 100.0 * (plscoeff[comp].VarX/tlv.VarX0);
    vxac += vx;
    sy = plscoeff[comp].SSYe[0]/tlv.SSY0[0];
  
    for (i=comp; i<=num_of_components; i++)
      r2[i] += (double) sy;
    
    sdec[comp] = sqrt((tlv.SSYp[0]-plscoeff[comp].SSYe[0])/x.nobj);
    sdec[comp] /= y.weights[0];
  }
  
  lsFreeTLV (&tlv);
 
  return 0;
}

void validate_pls (Matrix & x, Matrix & y, float * q2, float * sdep, 
    int numofcomp) 
{
  int i, comp;
  float ** ssy, ssy0;

  ssy = new float* [numofcomp+1];
  for (i=0; i<=numofcomp; i++) 
  {
    ssy[i] = new float;
    ssy[i][0] = 0.0;
  }
  
  ssy0 = 0.0;
    
  lsCheckMissing(&x,-99.00);
  lsCheckMissing(&y,-99.00);

  if ( (x.nmis + y.nmis) > 0 ) {
    lsSetMissing (&x, -99.00);
    lsSetMissing (&y, -99.00);
  }

  lsSetAverages (&y);

  for (i=0; i<y.nobj; i++) 
    ssy0 += sqr (y.data[i][0] - y.averages[0]);

  lsValidateLOO (&x,&y,numofcomp+1,1,ssy,
     NULL, NULL, NULL);

  /*
  lsValidateRG (&x, &y, numofcomp+1, 5, 20, 1, ssy, 
      NULL, NULL, NULL); 
  */

  for (comp=0; comp<=numofcomp; comp++) {
    q2[comp] = 1.00 - (ssy[comp][0]/ssy0);
    sdep[comp] = sqrt(ssy[comp][0]/(float)x.nobj);
  }

  for (i=0; i<=numofcomp; i++)
    delete ssy[i];
  delete[] ssy;

  return;
}

bool build_and_validate_pls (
    const char * prefix,
    std::vector<float> & logp_v, 
    std::vector<std::vector<int> > & finger_v, 
    int num_of_components,
    int max_num_of_components,
    bool validatepls)
{

  Matrix x, y;
  int numofobjs = fill_matrix (x, y, logp_v, finger_v);

  if (numofobjs > 10)
  {
    float * q2 = new float [max_num_of_components + 1];
    float * sdep = new float [max_num_of_components + 1];
    float * sdec = new float [max_num_of_components + 1];
    double * r2 = new double [max_num_of_components + 1];
    
    LV * plscoeff = new LV [max_num_of_components + 1];
                                                   
    for (int i=0; i<max_num_of_components + 1; i++) 
    {
      q2[i] = 0.0;
      sdep[i] = 0.0;
      sdec[i] = 0.0;
      r2[i] = 0.0;
    }
    
    build_pls (x, y, max_num_of_components, sdec, r2, plscoeff);

    fill_matrix (x, y, logp_v, finger_v);
    
    if (validatepls)
    {
      validate_pls (x, y, q2, sdep, max_num_of_components);

      std::cout << "n r2 q2 sdep" << std::endl;
      for (int j=0; j<=max_num_of_components; j++)
        std::cout << j+1 << " " << r2[j] << " " << 
          q2[j] << " " << sdep[j] << std::endl;
    }

    lsSetWeights (&x);
    lsSetWeights (&y);
                          
    // computing b0
    float zeroval = 0.0;
    std::vector<float> coeff;
    int max = num_of_components;
    for (int k=0; k<x.nvar; k++) 
    {
      coeff.push_back(plscoeff[max].b[0][k]*x.weights[k]/y.weights[0]);
      zeroval -= plscoeff[max].b[0][k]*x.weights[k]*x.averages[k]; 
    }
    zeroval /= y.weights[0];
    zeroval += y.averages[0];
    
    int counter = 0;
    std::stringstream ssp;
    ssp << prefix << "plscoeff.txt";
    std::cout << "Print Coeff " << ssp.str() << "... " << std::endl;
    if (access( ssp.str().c_str(), F_OK ) != -1)
      std::remove(ssp.str().c_str());
    std::ofstream fout (ssp.str().c_str(), std::ios::app);
   
    fout << "struct model2_pls prv_logp_"<<prefix<<"model = {" << std::endl;
    fout << std::scientific;
    fout << zeroval << ", { " << std::endl;
    std::vector<float>::iterator it = coeff.begin();
    for (; it != coeff.end(); ++it)
    {
      ++counter;
      fout << *it << ", ";
      if (counter%5 == 0)
        fout << std::endl;
    }
    fout << "}}; "<< std::endl;
    
    fout.close();
    
    lsFreeMatrix (&x); 
    lsFreeMatrix (&y);
    
    for (int j=0; j<max_num_of_components + 1; j++) 
      lsFreeLV (&plscoeff[j]);
    
    delete [] plscoeff;
                            
    delete [] q2;
    delete [] sdep;
    delete [] sdec;
    delete [] r2;
  }

  return true;
}
