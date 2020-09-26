//
//  main.cpp
//  Assignment6
//
//  Created by Minul Haque on 5/8/20.
//  Copyright © 2020 Minul Haque. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

// Initiallization routine for random number
void   RandomInitialise(int,int);
// Call to generate a random number between 0  and 1
double RandomUniform(void);
// Call to generate a Gaussian random number
double RandomGaussian(double,double);
// Call to generate a random integer between the two arguments
int    RandomInt(int,int);
// call to generate a random double between the two arguments
double RandomDouble(double,double);

int main() {
    
    RandomInitialise(1807,9373);
    int N = 128;
    
    double J = 0.1;
    while (J <= 0.85) {
        // Generate 2D lattice N x N
        int lattice[N][N];
        for (int i = 0; i < N; i ++) {
            for(int j = 0; j < N; j++) {
                
                // Generate spins at random with 0.65 positive bias
                 if (RandomUniform() < 0.65) {
                     lattice[i][j] = 1;
                 } else {
                     lattice[i][j] = -1;
                 }
            }
        }
        
        
        double preAvgE = 0;
        double preAvgEsq = 0;
        //double AvgE = 0;
        
        int i, j, k, l;
        double Si1, Si2;
        double Sj1, Sj2, Sl1, Sr1, Stop1, Sbot1, Sl2, Sr2, Stop2, Sbot2;
        double Hi1, Hi2, deltaH;

        // # MC Steps of N*N Moves
        int MCSteps = 1;
        while (MCSteps <= 100000) {
            double HiSum = 0.0;
            int moves = 1;
            while (moves <= N*N) {
                
                // Choose two random sites
                i = RandomInt(0, N-1);
                j = RandomInt(0, N-1);
                k = RandomInt(0, N-1);
                l = RandomInt(0, N-1);
                
                Si1 = lattice[i][j];
                Si2 = lattice[k][l];
                
                // Periodic Boundary Conditions for both sites
                Sl1 = lattice[i][(j - 1 + N) % N];
                Sr1 = lattice[i][(j + 1 + N) % N];
                Stop1 = lattice[(i - 1 + N) % N][j];
                Sbot1 = lattice[(i + 1 + N) % N][j];
                
                Sl2 = lattice[k][(l - 1 + N) % N];
                Sr2 = lattice[k][(l + 1 + N) % N];
                Stop2 = lattice[(k - 1 + N) % N][l];
                Sbot2 = lattice[(k + 1 + N) % N][l];
                
                Sj1 = Sl1 + Sr1 + Stop1 + Sbot1;
                Sj2 = Sl2 + Sr2 + Stop2 + Sbot2;
                
                Hi1 = -1.0 * J * Si1 * Sj1;
                Hi2 = -1.0 * J * Si2 * Sj2;
                
                deltaH = 2 * J * (Si1 * Sj1 + Si2 * Sj2);
                
                // Attempt Exchange if different
                if (lattice[i][j] != lattice[k][l]) {
                    if (deltaH < 0) {
                        lattice[i][j] = Si2;
                        lattice[k][l] = Si1;
                    } else if (RandomUniform() < exp(-1*deltaH)) {
                        lattice[i][j] = Si2;
                        lattice[k][l] = Si1;
                    }
                }
                moves++;
            }
            
            // Sum the energies of each site per MCStep
            if (MCSteps >= 20000) {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                                                                                    
                    Si1 = lattice[i][j];
                    Sl1 = lattice[i][(j - 1 + N) % N];
                    Sr1 = lattice[i][(j + 1 + N) % N];
                    Stop1 = lattice[(i - 1 + N) % N][j];
                    Sbot1 = lattice[(i + 1 + N) % N][j];
                                                                                    
                    Sj1 = Sl1 + Sr1 + Stop1 + Sbot1;

                    Hi1 = -1.0 * J * Si1 * Sj1;
                                                                                    
                    HiSum += Hi1;
                }
            }
            
            cout << MCSteps << "\t" << HiSum << endl;
            
            preAvgE += HiSum;
            preAvgEsq += pow(HiSum, 2);
            }
            MCSteps++;
        }
        
        double AvgE = preAvgE / 80000; // average energy
        double avgEsq = pow(AvgE, 2);
        double EsqAvg = preAvgEsq / 80000;
        double Cv = (EsqAvg - avgEsq);
       // cout << J << "\t" << Cv << endl;
        
        /*
        // Final Morphology for each J
         for (int i = 0; i < N; i++) {
             for (int j = 0; j < N; j++) {
                 cout << lattice[i][j] << "\t" << flush;
             }
             cout << endl;
         }
         cout << endl;
         cout << endl; */
        J+=0.1;
    }
}


#define FALSE 0
#define TRUE 1

/*
 
 This Random Number Generator is based on the algorithm in a FORTRAN
 version published by George Marsaglia and Arif Zaman, Florida State
 
 University; ref.: see original comments below.
 
 At the fhw (Fachhochschule Wiesbaden, W.Germany), Dept. of Computer
 Science, we have written sources in further languages (C, Modula-2
 Turbo-Pascal(3.0, 5.0), Basic and Ada) to get exactly the same test
 results compared with the original FORTRAN version.
 April 1989
 Karl-L. Noell <NOELL@DWIFH1.BITNET>
 and  Helmut  Weber <WEBER@DWIFH1.BITNET>
 
 This random number generator originally appeared in "Toward a Universal
 Random Number Generator" by George Marsaglia and Arif Zaman.
 Florida State University Report: FSU-SCRI-87-50 (1987)
 It was later modified by F. James and published in "A Review of Pseudo-
 random Number Generators"
 
 THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
 
 (However, a newly discovered technique can yield
 a period of 64^600. But that is still in the development stage.)
 It passes ALL of the tests for random number generators and has a period
 of 2^144, is completely portable (gives bit identical results on all
 machines with at least 24-bit mantissas in the floating point
 representation).
 
 The algorithm is a combination of a Fibonacci sequence (with lags of 97
 and 33, and operation "subtraction plus one, modulo one") and an
 "arithmetic sequence" (using subtraction).
 
 
 Use IJ = 1802 & KL = 9373 to test the random number generator. The
 subroutine RANMAR should be used to generate 20000 random numbers.
 Then display the next six random numbers generated multiplied by 4096*4096
 If the random number generator is working properly, the random numbers
 
 should be:
 
 6533892.0  14220222.0  7275067.0
 
 6172232.0  8354498.0   64633180.0
 
 */


/* Globals */

double u[97],c,cd,cm;

int i97,j97;

int test = FALSE;


/*
 
 This is the initialization routine for the random number generator.
 NOTE: The seed variables can have values between:    0 <= IJ <= 31328
 0 <= KL <= 30081
 
 The random number sequences created by these two seeds are of sufficient
 length to complete an entire calculation with. For example, if sveral
 different groups are working on different parts of the same calculation,
 each group could be assigned its own IJ seed. This would leave each group
 with 30000 choices for the second seed. That is to say, this random
 number generator can create 900 million different subsequences -- with
 each subsequence having a length of approximately 64^30.
 
 */

void RandomInitialise(int ij,int kl)

{
    
    double s,t;
    
    int ii,i,j,k,l,jj,m;
    
    
    /*
     
     Handle the seed range errors
     
     First random number seed must be between 0 and 31328
     
     Second seed must have a value between 0 and 30081
     
     */
    
    if (ij < 0 || ij > 31328 || kl < 0 || kl > 30081) {
        
        ij = 1802;
        
        kl = 9373;
        
    }
    
    
    i = (ij / 177) % 177 + 2;
    
    j = (ij % 177)       + 2;
    
    k = (kl / 169) % 178 + 1;
    
    l = (kl % 169);
    
    
    for (ii=0; ii<97; ii++) {
        
        s = 0.0;
        
        t = 0.5;
        
        for (jj=0; jj<24; jj++) {
            
            m = (((i * j) % 179) * k) % 179;
            
            i = j;
            
            j = k;
            
            k = m;
            
            l = (53 * l + 1) % 169;
            
            if (((l * m % 64)) >= 32)
                
                s += t;
            
            t *= 0.5;
            
        }
        
        u[ii] = s;
        
    }
    
    
    c    = 362436.0 / 16777216.0;
    
    cd   = 7654321.0 / 16777216.0;
    
    cm   = 16777213.0 / 16777216.0;
    
    i97  = 97;
    
    j97  = 33;
    
    test = TRUE;
    
}


/*
 
 This is the random number generator proposed by George Marsaglia in
 
 Florida State University Report: FSU-SCRI-87-50
 
 */

double RandomUniform(void)

{
    
    double uni;
    
    
    /* Make sure the initialisation routine has been called */
    
    if (!test)
        
        RandomInitialise(1802,9373);
    
    
    uni = u[i97-1] - u[j97-1];
    
    if (uni <= 0.0)
        
        uni++;
    
    u[i97-1] = uni;
    
    i97--;
    
    if (i97 == 0)
        
        i97 = 97;
    
    j97--;
    
    if (j97 == 0)
        
        j97 = 97;
    
    c -= cd;
    
    if (c < 0.0)
        
        c += cm;
    
    uni -= c;
    
    if (uni < 0.0)
        
        uni++;
    
    
    return(uni);
    
}


/*
 
 ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
 
 THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
 
 VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
 
 The function returns a normally distributed pseudo-random number
 
 with a given mean and standard devaiation.  Calls are made to a
 
 function subprogram which must return independent random
 
 numbers uniform in the interval (0,1).
 
 The algorithm uses the ratio of uniforms method of A.J. Kinderman
 
 and J.F. Monahan augmented with quadratic bounding curves.
 
 */

double RandomGaussian(double mean,double stddev)

{
    
    double  q,u,v,x,y;
    
    
    /*
     
     Generate P = (u,v) uniform in rect. enclosing acceptance region
     
     Make sure that any random numbers <= 0 are rejected, since
     
     gaussian() requires uniforms > 0, but RandomUniform() delivers >= 0.
     
     */
    
    do {
        
        u = RandomUniform();
        
        v = RandomUniform();
        
        if (u <= 0.0 || v <= 0.0) {
            
            u = 1.0;
            
            v = 1.0;
            
        }
        
        v = 1.7156 * (v - 0.5);
        
        
        /*  Evaluate the quadratic form */
        
        x = u - 0.449871;
        
        y = fabs(v) + 0.386595;
        
        q = x * x + y * (0.19600 * y - 0.25472 * x);
        
        
        /* Accept P if inside inner ellipse */
        
        if (q < 0.27597)
            
            break;
        
        
        /*  Reject P if outside outer ellipse, or outside acceptance region */
        
    } while ((q > 0.27846) || (v * v > -4.0 * log(u) * u * u));
    
    
    /*  Return ratio of P's coordinates as the normal deviate */
    
    return (mean + stddev * v / u);
    
}


/*
 
 Return random integer within a range, lower -> upper INCLUSIVE
 
 */

int RandomInt(int lower,int upper)

{
    
    return((int)(RandomUniform() * (upper - lower + 1)) + lower);
    
}


/*
 
 Return random float within a range, lower -> upper
 
 */

double RandomDouble(double lower,double upper)

{
    
    return((upper - lower) * RandomUniform() + lower);
    
}
