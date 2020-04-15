#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>

#include <cmath>

#include <getopt.h>
#include <unistd.h>

#include <vector>
#include <map>
#include <set>

#include <cassert>

#include <liblogd/logdpmanager.hpp>

#include <../src/logdp_private.h>

#include <../../grid_ng/inc/grinpp.h>
#include <../../grid_ng/inc/grinpp_grub.h>

#define HETDIM 146
#define HDVECDIM 118
extern grub_s_header grub_s_header_vec_static[HDVECDIM];

#ifdef USE_ALSO_FLAP
#include <../../flap-db/src/MFmolecule.h>
#endif 

#include <AtomAnalysis.h>

#include "common.hpp"

#define MAXNAMEDIM 1024

static inline void loadbar(int x, int n, int r, int w)
{
    // Only update r times.
    if ( x % (n/r +1) != 0 ) return;
 
    // Calculuate the ratio of complete-to-incomplete.
    float ratio = (float)x/(float)n;
    int   c     = (int) (ratio * w);
 
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    for (int x=0; x<c; x++)
      printf("=");
 
    for (int x=c; x<w; x++)
      printf(" ");
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
    printf("]\n\033[F\033[J");
}

void usagesex (char * name) 
{
  std::cerr << "usage: " << name << " [option] filenamei.mol2" << std::endl;
  std::cerr << " -h, --help                   : display this help and exit" << std::endl;
  std::cerr << "                                NOCTANOL = 1 CYCLEHEXANE = 2" << std::endl;
                           
  exit (1);
}

int
main (int argc, char **argv)
{
  int solvnum = 1;

  while (1) 
  {
    std::string inputs;
    std::vector<std::string>::iterator iter;
    std::vector<std::string> pairof_model_num, model_num;

    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"solvent", 0, NULL, 's'},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "hns:", long_options, &option_index);
    
    if (c == -1)
      break;

    switch (c) {
      case 's':
        solvnum = atoi(optarg);
        if (!((solvnum == 1) || (solvnum == 2)))
          usagesex (argv[0]);
        break;
      case 'h':
        usagesex (argv[0]);
        break;
      default:
        break;
    }
  }

  if (optind >= argc) 
    usagesex (argv[0]);

  std::string solvname = "";
  if (solvnum == 1)
    solvname = "noctanol";
  else if (solvnum == 2)
    solvname = "cyclehexane";

  std::set<std::string> ttypeset;
#ifdef USE_ALSO_FLAP
  std::set<std::string> flap_typeset;
#endif 

  struct logdp_moltype logdp_molecule;
  AA_Context * aacontext = aa3_create_context();

  char molname[MAXNAMEDIM];
  strncpy (molname, argv[optind], MAXNAMEDIM);

  int natoms = aa3_readFile (aacontext, ATOMANALYSIS_MOL2, molname, 
        MULTI_MOL2_NO, NULL);
  if (natoms == 0)
  {
    aa3_free (aacontext);
    aa3_destroy_context(aacontext);
    return 0;
  }

  aa3_setErrorLevel (aacontext, -1);
  aa3_dearomatize (aacontext);
  aa3_set_Noxide_bondtype (aacontext, 1);  

  if (logdp_aa_to_moltype (aacontext, &logdp_molecule, 0) == 0)
  {
    std::cout << molname << std::endl;
    for (int i = 0; i < logdp_molecule.dim; i++)
    { 
      if (strcmp (logdp_molecule.atom[i].ttype, "") != 0)
        ttypeset.insert(logdp_molecule.atom[i].ttype);
      fprintf (stdout, "%s \n", logdp_molecule.atom[i].ttype);
    }
  }

#ifdef USE_ALSO_FLAP
  struct pka_mols pka_molecola_in, pka_molecola;;
  if (pka_aa_to_mols (&pka_molecola_in, aacontext) == 0)
  {
    if (pka_move_hyd_to_end (&pka_molecola_in, &pka_molecola) == 0)
    {
      std::string * s_desc = NULL;
      // There could be some problem with H not at the end 
      // check with flap
      int res = get_atom_s_descriptors (&pka_molecola, s_desc, 1);
      int j = 0;
      for (int i = 0; i < pka_molecola.numofatom; i++)
      {
        if (strcmp (pka_molecola.atom[i], "H") != 0)
        {
          if (s_desc[j] != "")
            flap_typeset.insert(s_desc[j]);
          //std::cout << s_desc[j] << std::endl;
          ++j;
        }
      }
      delete [] s_desc;
    }
  }

  pka_free_mol (&pka_molecola);
  pka_free_mol (&pka_molecola_in);
#endif

  aa3_free (aacontext);
  aa3_destroy_context(aacontext);
  logdp_free_moltype (&logdp_molecule);

  return 0;
}
