/*
 * A model of protective mutualism
 * Uses a lattice with different productivity on each patch
 * Timothée Poisot 2013
 * t.poisot@gmail.com
 * Released under the terms of the GPL
 * Uses the GNU Scientific Library
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <getopt.h>
#include <unistd.h>

string FtS(double x)
{
   stringstream s;
   s << x;
   return s.str();
}

class patch
{
   public:
      // Population sizes
      double H;
      double P;
      double M;
      // Primary productivity
      double r;
      // Incoming
      double Hin;
      double Pin;
      double Min;
      // Outgoing
      double Hout;
      double Mout;
      double Pout;
};

int Wi = 80;
int He = 80;

double AverageR = 1.70;
double VarianceR = 1.35;

double HostDispersal = 0.01;
double SymbiontDispersal = 0.01;
double EnemyDispersal = 0.01;

int SimSteps = 5000;
double scalar = 0.005;
int OutSteps = 5;

// Parameters
double q = 0.005;
double B = 0.1;
double u = 1.9;
double a = 0.5;
double g = 0.1;
double de = 0.018;
double dm = 0.1;

int main(int argc, char *argv[])
{
   // Options
   int opt = 0;
   static struct option long_options[] = {
      {"rmean",   optional_argument, 0, 'r'},
      {"rvar",    optional_argument, 0, 'v'},
      {"hdisp",   optional_argument, 0, 'H'},
      {"pdisp",   optional_argument, 0, 'P'},
      {"mdisp",   optional_argument, 0, 'M'},
      {"alpha",   optional_argument, 0, 'a'},
      {0, 0, 0, 0}
   };
   int long_index = 0;
   while ((opt = getopt_long(argc, argv, "", long_options, &long_index))!=-1)
   {
      switch (opt)
      {
         case 'r':
            AverageR = atof(optarg);
            break;
         case 'v':
            VarianceR = atof(optarg);
            break;
         case 'H':
            HostDispersal = atof(optarg);
            break;
         case 'P':
            EnemyDispersal = atof(optarg);
            break;
         case 'M':
            SymbiontDispersal = atof(optarg);
            break;
         case 'a':
            a = atof(optarg);
            break;
         default:
            exit(EXIT_FAILURE);
      }
   }
   // Current time
   time_t begin_time = time(0);
   // RNG
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
   gsl_rng_set(rng, begin_time);
   // Out file
   string fname = "out-v"+FtS(VarianceR)+"-r"+FtS(AverageR)+"-H"+FtS(HostDispersal)+"-M"+FtS(SymbiontDispersal)+"-P"+FtS(EnemyDispersal)+"-a"+FtS(a)+"-i"+FtS(gsl_rng_uniform_pos(rng))+".dat";
   ofstream outfile;
   outfile.open(fname);
   outfile << "t x y r h p m\n";
   // The world
   patch W[Wi][He];
   // Set the world
   for (int x = 0; x < Wi ; ++x)
   {
      for(int y = 0; y < He; ++y)
      {
         W[x][y].r = AverageR + gsl_ran_gaussian(rng, VarianceR);
         W[x][y].r = W[x][y].r < 0.0 ? 0.0 : W[x][y].r;
         W[x][y].H = 10.0 + gsl_ran_gaussian(rng, 1.0);
         W[x][y].P = 1.0 + gsl_ran_gaussian(rng, 0.8);
         W[x][y].P = W[x][y].P < 0.0 ? 0.0 : W[x][y].P;
         W[x][y].M = 1.0 + gsl_ran_gaussian(rng, 0.8);
         W[x][y].M = W[x][y].M < 0.0 ? 0.0 : W[x][y].M;
         W[x][y].Hin = 0.0;
         W[x][y].Min = 0.0;
         W[x][y].Pin = 0.0;
         W[x][y].Hout = 0.0;
         W[x][y].Mout = 0.0;
         W[x][y].Pout = 0.0;
      }
   }
   // Start the iterations
   for(int ti = 0; ti <= SimSteps; ++ti)
   {
      // Population dynamics
      for(int x = 0; x < Wi; ++x)
      {
         for(int y = 0; y < He; ++y)
         {
            // Derivatives
            double dH = 0.0;
            double dP = 0.0;
            double dM = 0.0;
            // Temp variables
            double MutImpact = u/(u+W[x][y].M);
            // Dynamics
            dH += W[x][y].r - q*W[x][y].H - B*(W[x][y].P*MutImpact+a*W[x][y].M);
            dP += B*g*W[x][y].H*MutImpact - de;
            dM += B*g*a*W[x][y].H - dm;
            // Scaling
            dH *= W[x][y].H;
            dP *= W[x][y].P;
            dM *= W[x][y].M;
            // Update
            W[x][y].H += dH*scalar;
            W[x][y].P += dP*scalar;
            W[x][y].M += dM*scalar;
         }
      }
      // Migration
      for(int x = 0; x < Wi; ++x)
      {
         for(int y = 0; y < He; ++y)
         {
            int x_co[3];
            int y_co[3];
            x_co[0] = x-1;
            x_co[1] = x;
            x_co[2] = x+1;
            y_co[0] = y-1;
            y_co[1] = y;
            y_co[2] = y+1;
            for(int co = 0; co < 3; ++co)
            {
               x_co[co] = x_co[co] < 0 ? Wi + x_co[co] : x_co[co];
               y_co[co] = y_co[co] < 0 ? He + y_co[co] : y_co[co];
               x_co[co] = x_co[co] >= Wi ? x_co[co] - Wi : x_co[co];
               y_co[co] = y_co[co] >= He ? y_co[co] - He : y_co[co];
            }
            for(int xco = 0; xco < 3; ++xco)
            {
               int X = x_co[xco];
               for(int yco = 0; yco < 3; ++yco)
               {
                  int Y = y_co[yco];
                  if((x!=X)||(y!=Y))
                  {
                     // What goes out of the patch
                     W[x][y].Hout += W[x][y].H * scalar * HostDispersal * (1/8.0);
                     W[x][y].Pout += W[x][y].P * scalar * EnemyDispersal * (1/8.0);
                     W[x][y].Mout += W[x][y].M * scalar * SymbiontDispersal * (1/8.0);
                     // What goes into the patch
                     W[x][y].Hin += W[X][X].H * scalar * HostDispersal * (1/8.0);
                     W[x][y].Pin += W[X][Y].P * scalar * EnemyDispersal * (1/8.0);
                     W[x][y].Min += W[X][Y].M * scalar * SymbiontDispersal * (1/8.0);
                  }
               }
            }
         }
      }
      // Updates
      for(int x = 0; x < Wi; ++x)
      {
         for(int y = 0; y < He; ++y)
         {
            W[x][y].H += (W[x][y].Hin - W[x][y].Hout);
            W[x][y].P += (W[x][y].Pin - W[x][y].Pout);
            W[x][y].M += (W[x][y].Min - W[x][y].Mout);
            W[x][y].Hin = 0.0;
            W[x][y].Pin = 0.0;
            W[x][y].Min = 0.0;
            W[x][y].Hout = 0.0;
            W[x][y].Pout = 0.0;
            W[x][y].Mout = 0.0;
            // Print if needed
            if(ti % OutSteps == 0)
            {
               outfile << ti << " " << x << " " << y <<  " " << W[x][y].r  << " " << W[x][y].H << " " << W[x][y].P << " " << W[x][y].M << "\n";
            }
         }
      }
   }
   cout << "Execution complete in " << time(0) - begin_time << " seconds\n";
   outfile.close();
   gsl_rng_free(rng);
   return EXIT_SUCCESS;
}
