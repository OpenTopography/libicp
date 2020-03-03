/*
Copyright 2011. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

This file is part of libicp.
Authors: Andreas Geiger
Updated by Chelsea Scott

libicp is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or any later version.

libicp is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
libicp; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301, USA
*/

// Demo program showing how libicp can be used

#include <iostream>
#include "icpPointToPlane.h"
#include <list>
#include <array>
#include <fstream>  // std::ifstream
#include <iostream> // std::cout
#include <stdio.h> // printf/scanf family
#include "lasreader.hpp"
#include <vector>
#include <string>

using namespace std;

int main (int argc, char** argv) {

    if (argc != 3)  {
        printf("enter 3 arguments only eg. icp input_filename output_filename\n");
        return 0;
    }

    char* inputFileName = argv[1];
    char* outputFileName = argv[2];

    //output file is disp.txt
    //FILE* outfile = fopen("disp.txt","w");
    FILE* outfile = fopen(outputFileName,"w");

    //tiles.txt contains the coordinates of the core points
    //FILE* myfile = fopen("tiles.txt","r");
    FILE* myfile = fopen(inputFileName,"r");
    int x, y;

    //loop through tile.txt
    while (fscanf(myfile, "%d %d", &x,&y)==2){

        // write name of compare and reference tile to read
        char compare_file[50];
        char reference_file[50];
        snprintf(compare_file, 50,"las_diff/compare_%d_%d.las",x,y);
        snprintf(reference_file,50,"las_diff/reference_%d_%d.las",x,y);

        //read the compare tile
        LASreadOpener lasreadopener;
        lasreadopener.set_file_name(compare_file);
        LASreader* lasreader=lasreadopener.open();

        int npts = lasreader->npoints;
        double* M = (double*)malloc(npts*3*sizeof(double));
        int numM = npts;

        LASpoint* p;
        int i=0;

        while (lasreader -> read_point()) {
            p = &lasreader->point;
            M[i] = p->get_x();
            M[i+1] = p->get_y();
            M[i+2] = p->get_z();
            i += 3;
            }

        lasreader->close();
        delete lasreader;

        //read the reference file
        LASreadOpener nreadopener;
        nreadopener.set_file_name(reference_file);
        LASreader* nreader = nreadopener.open();
        npts=nreader->npoints;
        int numN = npts;

        double* N = (double*)malloc(npts*3*sizeof(double));
        i=0;
        while (nreader->read_point()){
            p= &nreader->point;
            N[i] = p-> get_x();
            N[i+1] = p->get_y();
            N[i+2] = p->get_z();
            i += 3;
        }

        nreader->close();
          delete nreader;

        // find the mean (x,y,z) of the compare dataset.
        double sum1=0;
        double sum2=0;
        double sum3=0;

        for (int k=0; k<numM; k+=1) {
            sum1 += M[k*3];
            sum2 += M[k*3+1];
            sum3 += M[k*3+2];
        }

        sum1 = sum1/numM;
        sum2 = sum2/numM;
        sum3 = sum3/numM;

        // subtract the mean from the reference and the compare
        for (int k=0; k<numM; k+=1) {
            M[k*3] = M[k*3] - sum1;
            M[k*3+1] = M[k*3+1] - sum2;
            M[k*3+2] = M[k*3+2] - sum3;
        }

        for (int k=0; k<numN; k+=1) {
            N[k*3] = N[k*3] - sum1;
            N[k*3+1] = N[k*3+1] - sum2;
            N[k*3+2] = N[k*3+2] - sum3;
        }

        // start with identity as initial transformation
        Matrix R = Matrix::eye(3);
        Matrix t(3,1);


        // run point-to-plane ICP (-1 = no outlier threshold)
        // Modified by Chelsea to add the "if" condition
        // orig: demo-bk-2020-01-21.cpp

        if (numN > 50 && numM > 50) {
            IcpPointToPlane icp(N,numN,3);
            double residual = icp.fit(M,numM,R,t,-1);

            // write the output file
            fprintf(outfile, "%d %d %f %f %f %f %f %f %f %f %f %f %f %f\n",x,y, t.val[0][0], t.val[1][0], t.val[2][0], R.val[0][0], R.val[0][1], R.val[0][2], R.val[1][0], R.val[1][1], R.val[1][2], R.val[2][0], R.val[2][1], R.val[2][2] );
        }
        // free memory
        free(M);
        free(N);
    }
    fclose(outfile);
    fclose(myfile);
    // success
    return 0;
}
