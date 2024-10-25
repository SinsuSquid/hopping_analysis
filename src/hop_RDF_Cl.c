#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMBINS 1000
#define NUMINTERVALS 6
#define DELTA_T 1

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;
double box_x, box_y, box_z;

int numLi, numCl, numAl;

int readTraj(void);

FILE *fp_in;
FILE *fp_out;

// distList[atomIdx][start][t]
int deltaT;
double ***avgPosA;
double ***avgPosB;
double **hop;
double getDistance(double *, double *);
void initialize();
void hop_function();
void hop_RDF();

int main(int argc, char *argv[]){
	if (argc != 4){
		printf("USAGE : ./hop.x ***.lammpstrj DELTA_T hop_RDF.out\n");
		exit(1);
	}
	fp_in = fopen(argv[1], "r");
	deltaT = atoi(argv[2]);
	fp_out = fopen(argv[3], "w");

	atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
	coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

	readTraj();

	printf("\n\t\t< FILE SUMMARY >\n\n");
	printf("\tNumber of Timesteps : %d\n", numTraj);
	printf("\tNumber of Atoms : %d\n\n", numAtoms[0]);

	// Don't use this w/ NpT !
	box_x = (box[0][0][1] - box[0][0][0]);
	box_y = (box[0][1][1] - box[0][1][0]);
	box_z = (box[0][2][1] - box[0][2][0]);

	numLi = 0; numCl = 0; numAl = 0;
	for (int i = 0; i < numAtoms[0]; i++){
		switch(atom[0][i][1]){
			case 1: numLi += 1; break;
			case 2: numCl += 1; break;
			case 3: numAl += 1; break;
		}
	}
	printf("\tbox_x = %lf\n", box_x * 2.0);
	printf("\tbox_y = %lf\n", box_y * 2.0);
	printf("\tbox_z = %lf\n", box_z * 2.0);
	printf("\tnumLi = %d\n", numLi);
	printf("\tnumCl = %d\n", numCl);
	printf("\tnumAl = %d\n", numAl);
	printf("\n");

	initialize();
	hop_function();
	hop_RDF();

	free(atom);
	free(coord);
	free(hop);

	printf("\tAll Tasks are Done ! >:D\n");

	return 0;
}

void hop_RDF(){
	printf("\tCalculating RDF ...\n");
	fprintf(fp_out, "#\tr\ttotal\t0.0-1.0-\t1.0-2.0\t2.0-3.0\t3.0-4.0\t4.0-5.0\n");

	double **RDF = (double **)malloc(sizeof(double *) * NUMINTERVALS);
	long *numPoints = (long *)malloc(sizeof(long) * NUMINTERVALS);
	for (int i = 0; i < NUMINTERVALS; i++){
		RDF[i] = (double *)malloc(sizeof(double) * NUMBINS);
		numPoints[i] = (long)0;
		for (int ii = 0; ii < NUMBINS; ii++) RDF[i][ii] = 0.0;
	}
	double binsize = sqrt((box_x * box_x) + (box_y * box_y) + (box_z * box_z)) / NUMBINS;

	for (int t = deltaT; t < numTraj - deltaT; t++){
		if (!(t % 100)) printf("\t\tNow t = %d...\n", t);
		for (int i = 0; i < numLi; i++){
			int hopIdx = round(hop[t][i]);
			for (int j = numLi; j < numLi + numCl; j++){
				double distance = sqrt(getDistance(coord[t][i], coord[t][j]));
				int idx = (int)(distance / binsize);
				if (idx < NUMBINS){
					RDF[0][idx] += 2; numPoints[0] += 2;
					if (hopIdx <= 4){
						RDF[1 + hopIdx][idx] += 2; numPoints[1 + hopIdx] += 2;
					}
				}
			}
		}
	}

	double volume = box_x * box_y * box_z;
	for (int i = 0; i < NUMBINS; i++){
		double r = i * binsize;
		fprintf(fp_out, "%lf", r);
		for (int ii = 0; ii < NUMINTERVALS; ii++){
			double val = RDF[ii][i];
			val /= (numPoints[ii] / volume);
			val /= 4.0 * M_PI * r * r * binsize;
			fprintf(fp_out, "\t%lf", val);
		}
		fprintf(fp_out, "\n");
	}
	
	return;
}

double getDistance(double *coord1, double *coord2){
	double dx = coord1[0] - coord2[0];
	double dy = coord1[1] - coord2[1];
	double dz = coord1[2] - coord2[2];

	dx -= round(dx / box_x) * box_x;
	dy -= round(dy / box_y) * box_y;
	dz -= round(dz / box_z) * box_z;

	return dx * dx + dy * dy + dz * dz;
}

void initialize(){
	printf("\tNow Initializing Averaged Position ...\n");

	int GAP = deltaT / 2; printf("\tGAP = %d\n", GAP);

	avgPosA = (double ***)malloc(sizeof(double **) * numTraj);
	avgPosB = (double ***)malloc(sizeof(double **) * numTraj);

	for (int t = GAP; t < numTraj - GAP; t++){
		double **tempPosA = (double **)malloc(sizeof(double *) * numLi);
		double **tempPosB = (double **)malloc(sizeof(double *) * numLi);

		for (int idx = 0; idx < numLi; idx++){
			double *pA = (double *)malloc(sizeof(double) * 3);
			double *pB = (double *)malloc(sizeof(double) * 3);

			for (int dt = 0; dt < GAP; dt++){
				pA[0] += coord[t - dt][idx][0]; pA[1] += coord[t - dt][idx][1]; pA[2] += coord[t - dt][idx][2];
				pB[0] += coord[t + dt][idx][0]; pB[1] += coord[t + dt][idx][1]; pB[2] += coord[t + dt][idx][2];
			}

			for (int i = 0; i < 3; i++) { pA[i] /= GAP; pB[i] /= GAP; }
			// printf("pA[%d][%d] = [%lf, %lf, %lf]\t pB[%d][%d] = [%lf, %lf, %lf]\n", t, idx, pA[0], pA[1], pA[2], t, idx, pB[0], pB[1], pB[2]);

			tempPosA[idx] = pA;
			tempPosB[idx] = pB;
		}
		avgPosA[t] = tempPosA;
		avgPosB[t] = tempPosB;
	} 
}

void hop_function(){
	printf("\tNow Calculating Hop Function ...\n");
	int GAP = deltaT / 2; printf("\tGAP = %d\n", GAP);
	// fprintf(fp_out, "#\tt\tatomIdx\n");

	hop = (double **)malloc(sizeof(double *) * numTraj);
	for (int t = GAP + GAP; t < numTraj - GAP - GAP; t++){
		double *perTime = (double *)malloc(sizeof(double) * numLi);

		box_x = (box[t][0][1] - box[t][0][0]);
		box_y = (box[t][1][1] - box[t][1][0]);
		box_z = (box[t][2][1] - box[t][2][0]);

		for (int idx = 0; idx < numLi; idx++){
			double dA = 0.0; double dB = 0.0;
			for (int dt = 0; dt < GAP; dt++){
				dB += getDistance(coord[t + dt][idx], avgPosA[t + dt][idx]);
				dA += getDistance(coord[t - dt][idx], avgPosB[t - dt][idx]);
			}
			dA /= GAP; dB /= GAP;
			// printf("dA[%d][%d] = %lf\tdB[%d][%d] = %lf\n", t, idx, dA, t, idx, dB);
			perTime[idx] = sqrt(dA * dB);
		}
		hop[t] = perTime;
		/*
		fprintf(fp_out, "%d\t", t);
		for (int i = 0; i < numLi; i++){
			fprintf(fp_out, "%lf\t",  hop[t][i]);
		}
		fprintf(fp_out, "\n");
		*/
	}
	printf("\tInitialized Hop Table !\n");
}

int readTraj(void){
	char *iostat;
	char line[LINESIZE];

	while(1){
		iostat = fgets(line, LINESIZE, fp_in);
		// printf("s", line);

		if (!iostat) break;
		else if (strcmp(line, "ITEM: TIMESTEP\n") == 0){
			int temp;
			fscanf(fp_in, "%d", &temp);
			timestep[numTraj] = temp;
			// printf("timestep : %d\n", timestep[numTraj]);
		}
		else if (strcmp(line, "ITEM: NUMBER OF ATOMS\n") == 0){
			int temp;
			fscanf(fp_in, "%d", &temp);
			numAtoms[numTraj] = temp;
			// printf("numAtoms : %d\n", numAtoms[numTraj]);
		}
		else if (strcmp(line, "ITEM: BOX BOUNDS pp pp pp\n") == 0){
			double temp1, temp2;
			for (int i = 0; i < 3; i++){
				fscanf(fp_in, "%lf %lf", &temp1, &temp2);
				box[numTraj][i][0] = temp1;
				box[numTraj][i][1] = temp2;
			}
			/*
			printf("box :\n%f %f\n%f %f\n%f %f\n",
			        box[numTraj][0][0], box[numTraj][0][1],
			        box[numTraj][1][0], box[numTraj][1][1],
			        box[numTraj][2][0], box[numTraj][2][1]);
			*/
		}
		else if (strcmp(line, "ITEM: ATOMS id type x y z\n") == 0){
			double x, y, z;
			int type, id, ix, iy, iz;

			double boxlength[3] = {box[numTraj][0][1] - box[numTraj][0][0],
					       box[numTraj][1][1] - box[numTraj][1][0],
					       box[numTraj][2][1] - box[numTraj][2][0]};

			int **atomPerTraj;
			atomPerTraj = (int**)malloc(sizeof(int*) * numAtoms[numTraj]);
			double **coordPerTraj;
			coordPerTraj = (double**)malloc(sizeof(double*) * numAtoms[numTraj]);

			for (int i = 0; i < numAtoms[numTraj]; i++){
				int *atomTemp;
				atomTemp = (int*)malloc(sizeof(int) * 2);
				double *coordTemp;
				coordTemp = (double*)malloc(sizeof(double) * 3);

				fscanf(fp_in, "%d %d %lf %lf %lf",
				       &id, &type, &x, &y, &z);

				atomTemp[0] = id;
				atomTemp[1] = type;
				atomPerTraj[i] = atomTemp;

				coordTemp[0] = x;
				coordTemp[1] = y;
				coordTemp[2] = z;
				coordPerTraj[i] = coordTemp;

				/*
				printf("timestep : %d id : %d, type : %d\n",
				       timestep[numTraj], atomPerTraj[i][0], atomPerTraj[i][1]);
				printf("\tx : %f y : %f z : %f\n",
				       coordPerTraj[i][0],
				       coordPerTraj[i][1],
				       coordPerTraj[i][2]);
				*/
			}
			atom[numTraj] = atomPerTraj;
			coord[numTraj] = coordPerTraj;
			numTraj += 1;
		}
	}
	fclose(fp_in);
	return 0;
}
