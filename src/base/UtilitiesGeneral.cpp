#include "headersBasic.h"
#include "UtilitiesGeneral.h"


void SetTimeParametersFluid(int tis, double rho, double dt, VectorXd& td)
{
  td.setZero();

  double  alpf, alpm, beta, gamm;

  td[0] = dt;

  switch(tis)
  {
      case  0: // quasi static

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

       break;

      case  1: // generalised midpoint rule

            alpf = 1.0/(1.0 + rho);
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = alpf*dt; //af*dt
            td[6]  = alpm/gamm;
            td[7]  = 1.0-alpm/gamm;
            td[8]  = 0.0;

            td[9]  = 1.0/dt;  // v_{n+1}
            td[10] = -td[9];  // v_n
            td[11] = 0.0;     // v_{n-1}
            td[12] = 0.0;     // v_{n-2}
            td[15] = 0.0;     // a_n

        break;

      case  2: // generalised alpha-method

            alpf = 1.0/(1.0 + rho);
            alpm = 0.5*(3.0 - rho)/(1.0 + rho);
            gamm = 0.5 + alpm - alpf;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = alpf*dt;
            td[6]  = alpm/gamm;
            td[7]  = 1.0-alpm/gamm;
            td[8]  = alpm/gamm/dt;

            td[9]  = 1.0/gamm/dt;     // v_{n+1}
            td[10] = -td[9];          // v_n
            td[11] = 0.0;     // v_{n-1}
            td[12] = 0.0;     // v_{n-2}
            td[15] = 1.0 - 1.0/gamm;  // a_n

         break;

      case  3: // Backward Euler or BDF1

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

            td[5]  = dt;
            td[6]  = 1.0;
            td[7]  = 0.0;
            td[8]  = 1.0/dt;

            td[9]  =  1.0/dt;    // v_{n+1}
            td[10] = -1.0/dt;    // v_n
            td[11] =  0.0;       // v_{n-1}
            td[12] =  0.0;       // v_{n-2}
            td[13] =  0.0;       // v_{n-3}
            td[15] =  0.0;       // a_n

         break;

      case  4: // generalised alpha-method

            alpf = 1.0/(1.0 + rho);
            alpm = 0.5*(3.0 - rho)/(1.0 + rho);
            gamm = 0.5 + alpm - alpf;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = alpf*dt;
            td[6]  = alpm/gamm;
            td[7]  = 1.0-alpm/gamm;
            td[8]  = alpm/gamm/dt;

            td[9]  = 1.0/gamm/dt;    // v_{n+1}
            td[10] = -td[9];         // v_n
            td[11] = 0.0;            // v_{n-1}
            td[12] = 0.0;            // v_{n-2}
            td[13] =  0.0;           // v_{n-3}
            td[15] = 1.0 - 1.0/gamm; // a_n

         break;

      case  5: // BDF2

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

            td[5]  = dt;
            td[6]  = 1.0;
            td[7]  = 0.0;
            td[8]  = 1.0/dt;

            td[9]  =  1.5/dt;    // v_{n+1}
            td[10] = -2.0/dt;    // v_n
            td[11] =  0.5/dt;    // v_{n-1}
            td[12] =  0.0;       // v_{n-2}
            td[13] =  0.0;       // v_{n-3}
            td[15] =  0.0;       // a_n

         break;

      case  6: // BDF3

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

            td[5]  = dt;
            td[6]  = 1.0;
            td[7]  = 0.0;
            td[8]  = 1.0/dt;

            td[9]  =  (11.0/6.0)/dt;    // v_{n+1}
            td[10] = -3.0/dt;           // v_n
            td[11] =  1.5/dt;           // v_{n-1}
            td[12] = -(1.0/3.0)/dt;     // v_{n-2}
            td[13] =  0.0;              // v_{n-3}
            td[15] =  0.0;              // a_n

         break;

      case  7: // BDF4

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

            td[5]  = dt;
            td[6]  = 1.0;
            td[7]  = 0.0;
            td[8]  = 1.0/dt;

            td[9]  =  (25.0/12.0)/dt;    // v_{n+1}
            td[10] = -4.0/dt;           // v_n
            td[11] =  3.0/dt;           // v_{n-1}
            td[12] = -(4.0/3.0)/dt;     // v_{n-2}
            td[13] =  (1.0/4.0)/dt;     // v_{n-3}
            td[15] =  0.0;              // a_n

         break;

      default:
            cerr << " SetTimeParametersFluid ... invalid value of tis!" << endl;

         break;
  }

  return;
}


