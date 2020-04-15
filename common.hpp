#include <basStat.h>

int fill_matrix (Matrix & x, Matrix & y, std::vector<float> & logp_v, 
    std::vector<std::vector<int> > & finger_v);

int build_pls (Matrix & x, Matrix & y, int num_of_components, 
    float * sdec, double * r2, LV * plscoeff);

void validate_pls (Matrix & x, Matrix & y, float * q2, float * sdep, 
    int numofcomp);

bool build_and_validate_pls (
    const char * prefix,
    std::vector<float> & logp_v, 
    std::vector<std::vector<int> > & finger_v, 
    int num_of_components,
    int max_num_of_components, 
    bool validatepls);
