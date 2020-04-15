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
  std::cerr << "usage: " << name << " [option] database.txt" << std::endl;
  std::cerr << " -h, --help                   : display this help and exit" << std::endl;
  std::cerr << " -s, --solvent=[solv]         : specify solvent name (default=1)" << std::endl;
  std::cerr << "                                NOCTANOL = 1 CYCLEHEXANE = 2" << std::endl;
  std::cerr << " -n, --no-validate            : do not validate PLS" << std::endl;
                           
  exit (1);
}

int
main (int argc, char **argv)
{
  int solvnum = 1;
  bool validatepls = true;

  while (1) 
  {
    std::string inputs;
    std::vector<std::string>::iterator iter;
    std::vector<std::string> pairof_model_num, model_num;

    int c, option_index;
    static struct option long_options[] = {
      {"help", 0, NULL, 'h'},
      {"solvent", 0, NULL, 's'},
      {"no-validate", 0, NULL, 'n'},
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
      case 'n':
        validatepls = false;
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

  std::ifstream fin (argv[optind]);
  std::string molname;
  float logpval = 0.0;

  int numofline = 0;
  while (!fin.eof()) 
  {
    fin.ignore (1024, '\n');
    numofline++;
  }
  numofline--;
  fin.clear();
  fin.seekg(0);

  std::cout << "Total amount of lines: " << numofline << std::endl; 

  std::set<std::string> ttypeset;
#ifdef USE_ALSO_FLAP
  std::set<std::string> flap_typeset;
#endif 

  int conta = 0;
  std::cout << "Start producing AT list ..." << std::endl;
  do 
  {
    fin >> molname;
    fin >> logpval;

    struct logdp_moltype logdp_molecule;
    AA_Context * aacontext = aa3_create_context();

    char mname[MAXNAMEDIM];
    strncpy (mname, molname.c_str(), MAXNAMEDIM);

    int natoms = aa3_readFile (aacontext, ATOMANALYSIS_MOL2, mname, 
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
      //std::cout << molname << " " << logpval << std::endl;
      for (int i = 0; i < logdp_molecule.dim; i++)
      { 
        if (strcmp (logdp_molecule.atom[i].ttype, "") != 0)
          ttypeset.insert(logdp_molecule.atom[i].ttype);
        //fprintf (stdout, "%s \n", logdp_molecule.atom[i].ttype);
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

    ++conta;
    loadbar(conta, numofline, numofline/10, 80);
  } while (conta < numofline);

  unsigned int ttypedim = 0;
  // copialo in un file di header da usare poi per la numerazione 
  std::cout << "Writing atomtype file..." << std::endl;
  std::ofstream outf;
  outf.open("liblogdp_atomtype.h");
  outf << "#ifndef PRV_LIBLOGDP_ATOMTYPE_INC" << std::endl;
  outf << "#define PRV_LIBLOGDP_ATOMTYPE_INC" << std::endl;
  outf << "#define LIBLOGDP_ATOMTYPE_DIM " << ttypeset.size() << std::endl;
  outf << std::endl;
  outf << "struct liblogdp_atomtype_str prv_liblogdp_atomtype_stri" << 
    solvname << "[LIBLOGDP_ATOMTYPE_DIM] = {" 
       << std::endl;
  std::set<std::string>::const_iterator ttypeit = ttypeset.begin();
  for (; ttypeit != ttypeset.end(); ++ttypeit)
  {
    if (ttypedim < ttypeit->size() )
      ttypedim = ttypeit->size();
    outf << "  {\"" << *ttypeit << "\"}," << std::endl;
  }
  outf << "};" << std::endl;
  outf << "#endif" << std::endl;
  outf.close();

  ++ttypedim;

  std::cout << "\e[1m new TTYPEDIM " << ttypedim << "\e[0m" << std::endl;

  fin.clear();
  fin.seekg(0);

#ifdef USE_ALSO_FLAP  
  unsigned int flap_typedim = 0;
  // copialo in un file di header da usare poi per la numerazione 
  std::cout << "Writing flap_atomtype file..." << std::endl;
  std::ofstream foutf;
  foutf.open("liblogdp_flap_atomtype.h");
  foutf << "#ifndef PRV_LIBLOGDP_FLAP_ATOMTYPE_INC" << std::endl;
  foutf << "#define PRV_LIBLOGDP_FLAP_ATOMTYPE_INC" << std::endl;
  foutf << "#define LIBLOGDP_FLAP_ATOMTYPE_DIM " << flap_typeset.size() << std::endl;
  foutf << std::endl;
  foutf << "struct liblogdp_atomtype_fstr prv_liblogdp_atomtype_fstri" << 
    solvname << "[LIBLOGDP_FLAP_ATOMTYPE_DIM] = {" 
       << std::endl;
  std::set<std::string>::const_iterator flap_typeit = flap_typeset.begin();
  for (; flap_typeit != flap_typeset.end(); ++flap_typeit)
  {
    if (flap_typedim < flap_typeit->size() )
      flap_typedim = flap_typeit->size();
    foutf << "  {\"" << *flap_typeit << "\"}," << std::endl;
  }
  foutf << "};" << std::endl;
  foutf << "#endif" << std::endl;
  foutf.close();

  ++flap_typedim;

  std::cout << "\e[1m new FLAP_TYPEDIM " << flap_typedim << "\e[0m" << std::endl;

  fin.clear();
  fin.seekg(0);
#endif 

  std::vector<float> logpvals;
  std::vector<std::string> molnames;
  std::vector<std::vector<int> > modelttype;
  std::vector<std::vector<int> > modelgtype;
#ifdef USE_ALSO_FLAP  
  std::vector<std::vector<int> > modelftype;
#endif

  conta = 0;
  std::cout << "Start producing desc ..." << std::endl;
  do 
  {
    fin >> molname;
    fin >> logpval;

    struct logdp_moltype logdp_molecule;
    AA_Context * aacontext = aa3_create_context();

    char mname[MAXNAMEDIM];
    strncpy (mname, molname.c_str(), MAXNAMEDIM);

    int natoms = aa3_readFile (aacontext, ATOMANALYSIS_MOL2, mname, 
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

    std::vector<int> desc;
    for (int i=0; i<(int)ttypeset.size(); ++i)
      desc.push_back(0);

    std::vector<int> descg;
    for (int i=0; i<(int)HETDIM; ++i)
      descg.push_back(0);

    if (logdp_aa_to_moltype (aacontext, &logdp_molecule, 0) == 0)
    {
      //std::cout << molname << " " << logpval << std::endl;
      for (int i = 0; i < logdp_molecule.dim; i++)
      { 
        if (strcmp (logdp_molecule.atom[i].ttype, "") != 0)
        {
          std::set<std::string>::iterator find_result = ttypeset.find(
            logdp_molecule.atom[i].ttype);
          int idx = std::distance(ttypeset.begin(), find_result);
          desc[idx] = desc[idx] + 1;
        }

        int gidx = get_num_gtype(logdp_molecule.atom[i].gtype);
        if (gidx >= 0)
          descg[gidx] = descg[gidx] + 1;
      }

      modelttype.push_back(desc);
      logpvals.push_back(logpval);      
      molnames.push_back(molname);
    }

#ifdef USE_ALSO_FLAP    
    std::vector<int> flap_desc;
    for (int i=0; i<(int)flap_typeset.size(); ++i)
      flap_desc.push_back(0);

    struct pka_mols pka_molecola, pka_molecola_in;
    if (pka_aa_to_mols (&pka_molecola_in, aacontext) == 0)
    {
      if (pka_move_hyd_to_end (&pka_molecola_in, &pka_molecola) == 0)
      {
        std::string * s_desc = NULL;
        int res = get_atom_s_descriptors (&pka_molecola, s_desc, 1);
        int j = 0;
        for (int i = 0; i < pka_molecola.numofatom; i++)
        {
          if (strcmp (pka_molecola.atom[i], "H") != 0)
          {
            if (s_desc[j] != "")
            {
              std::set<std::string>::iterator find_result = flap_typeset.find(
                s_desc[j]);
              int idx = std::distance(flap_typeset.begin(), find_result);
              flap_desc[idx] = flap_desc[idx] + 1;
            }
            //std::cout << s_desc[j] << std::endl;
            ++j;
          }
        }
        delete [] s_desc;
        
        modelftype.push_back(flap_desc);
      }
    }
    pka_free_mol (&pka_molecola);
    pka_free_mol (&pka_molecola_in);
#endif

    if (aa3_assignGridTypes(aacontext) == 0)
    {
      for (int i = 0; i < logdp_molecule.dim; i++)
      { 
        int gidx = get_num_gtype(logdp_molecule.atom[i].gtype);
        if (gidx >= 0)
          descg[gidx] = descg[gidx] + 1;
      }

      modelgtype.push_back(descg);
    }

    aa3_free (aacontext);
    aa3_destroy_context(aacontext);
    logdp_free_moltype (&logdp_molecule);

    ++conta;
    loadbar(conta, numofline, numofline/10, 80);
  } while (conta < numofline);

  assert(modelttype.size() == logpvals.size());
  assert(modelgtype.size() == logpvals.size());
  assert(modelttype.size() == modelgtype.size());
  assert(modelttype.size() == molnames.size());
#ifdef USE_ALSO_FLAP  
  assert(modelftype.size() == molnames.size());
  assert(modelftype.size() == modelgtype.size());
#endif

  std::cout << "Writing desc files..." << std::endl;
  std::ofstream gdescfile, tdescfile, fdescfile;
  gdescfile.open((solvname+"_modelgtype.txt").c_str()); 
  tdescfile.open((solvname+"_modelttype.txt").c_str()); 

#ifdef USE_ALSO_FLAP  
  fdescfile.open((solvname+"_modelftype.txt").c_str());
#endif

  for (int i=0; i<(int)logpvals.size(); ++i)
  {
    gdescfile << molnames[i] << " ";
    for (int j=0; j<(int)modelgtype[i].size(); ++j)
      gdescfile << modelgtype[i][j] << " ";
    gdescfile << logpvals[i] << std::endl;

    tdescfile << molnames[i] << " ";
    for (int j=0; j<(int)modelttype[i].size(); ++j)
      tdescfile << modelttype[i][j] << " ";
    tdescfile << logpvals[i] << std::endl;

#ifdef USE_ALSO_FLAP
    fdescfile << molnames[i] << " ";
    for (int j=0; j<(int)modelftype[i].size(); ++j)
      fdescfile << modelftype[i][j] << " ";
    fdescfile << logpvals[i] << std::endl;
#endif

  } 
  gdescfile.close(); 
  tdescfile.close(); 
#ifdef USE_ALSO_FLAP
  fdescfile.close();
#endif

  std::cout << "Build and validate modelgtype..." << std::endl;
  build_and_validate_pls ((solvname+"_gtype_").c_str(), 
      logpvals, modelgtype, 4, 6, validatepls);

  std::cout << "Build and validate modelttype..." << std::endl;
  build_and_validate_pls ((solvname+"_ttype_").c_str(), 
      logpvals, modelttype, 25, 27, validatepls);

#ifdef USE_ALSO_FLAP
  std::cout << "Build and validate modelftype..." << std::endl;
  build_and_validate_pls ((solvname+"_ftype_").c_str(), 
      logpvals, modelftype, 25, 27, validatepls);
#endif

  return 0;
}

/*
ttype gtype piu' tipi direttamente legati conto
 1 per ogni atomo per ogni molecola conti il numero di gtype 
 200 gtpe data la moloecola A 1 di gtype 4 e 5 di gtype 6  
   0 0 0 1 0 5 0 0 0 0   

due modelli 

intercetta zero 
*/
