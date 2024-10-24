#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BUFFERSIZE 1000

typedef struct{
	float x, y, z;
} dcdCoord;

typedef struct{
	int id, type;
} particleInfo;

FILE *fp_init, *fp_dcd, *fp_out;

void initReader();
void dcdHeader();
void dcdReader(dcdCoord *, double *);
void outputWriter(dcdCoord *, double *, int);
int compare (const void * a, const void * b);

int numAtom, numTraj, numType;
char *title;
double *masses;
particleInfo *init;

int main(int argc, char *argv[]){
	if (argc != 4){
		for (int i = 0; i < argc; i++)
			printf("argv[%d] : %s\n", i, argv[i]);
		printf("USAGE : ./dcd2lmptrj input.lammps_data trajectory.dcd output.lammpstrj\n");
		exit(1);
	}

	fp_init = fopen(argv[1], "r");
	fp_dcd = fopen(argv[2], "rb");
	fp_out = fopen(argv[3], "w");

	title = (char *)malloc(sizeof(char) * BUFFERSIZE);

	dcdHeader();
	initReader();

	printf("\tTitle : %s\n", title);

	for (int i = 0; i < numTraj; i++){
		if (!(i % 100))
			printf("Now %d (st/nd/th) trajectory ...\n", i);
		double *box = (double *)malloc(sizeof(double) * 6);
		dcdCoord *coord = (dcdCoord *)malloc(sizeof(*coord) * numAtom);
		dcdReader(coord, box);
		outputWriter(coord, box, i);
	}

	fclose(fp_init);
	fclose(fp_dcd);
	fclose(fp_out);

	printf("Done! >:D\n");

	return 0;
}

void outputWriter(dcdCoord *coord, double *box, int idx){
	fprintf(fp_out, "ITEM: TIMESTEP\n");
	fprintf(fp_out, "%d\n", idx);
	fprintf(fp_out, "ITEM: NUMBER OF ATOMS\n");
	fprintf(fp_out, "%d\n", numAtom);
	fprintf(fp_out, "ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(fp_out, "%lf %lf\n", -0.5 * box[0], 0.5 * box[0]);
	fprintf(fp_out, "%lf %lf\n", -0.5 * box[2], 0.5 * box[2]);
	fprintf(fp_out, "%lf %lf\n", -0.5 * box[5], 0.5 * box[5]);
	fprintf(fp_out, "ITEM: ATOMS id type x y z\n");
	for (int i = 0; i < numAtom; i++){
		fprintf(fp_out, "%d %d %f %f %f\n",
				init[i].id, init[i].type,
				coord[i].x, coord[i].y, coord[i].z);
	}

	return;
}

void dcdReader(dcdCoord *coord, double *box){
	int *dummyInt = (int*)malloc(sizeof(int) * BUFFERSIZE);
	double *dummyDouble = (double*)malloc(sizeof(double) * BUFFERSIZE);
	float *dummyFloat = (float*)malloc(sizeof(float) * BUFFERSIZE);

	int vectorByte = numAtom * sizeof(float);

	fread(dummyInt, sizeof(int), 1, fp_dcd); // 48
	// dim = {a, gamma, b, beta, alpha, c}

	for (int i = 0; i < 6; i++){
		fread(dummyDouble, sizeof(double), 1, fp_dcd); // box
		box[i] = *dummyDouble;
	}

	fread(dummyInt, sizeof(int), 1, fp_dcd); // 48
	fread(dummyInt, sizeof(int), 1, fp_dcd); // vectorByte

	for (int i = 0; i < numAtom; i++){ // x
		fread(dummyFloat, sizeof(float), 1, fp_dcd);
		coord[i].x = *dummyFloat;
	}

	fread(dummyInt, sizeof(int), 1, fp_dcd); // vectorByte
	fread(dummyInt, sizeof(int), 1, fp_dcd); // vectorByte
	
	for (int i = 0; i < numAtom; i++){ // y
		fread(dummyFloat, sizeof(float), 1, fp_dcd);
		coord[i].y = *dummyFloat;
	}

	fread(dummyInt, sizeof(int), 1, fp_dcd); // vectorByte
	fread(dummyInt, sizeof(int), 1, fp_dcd); // vectorByte

	for (int i = 0; i < numAtom; i++){ // z
		fread(dummyFloat, sizeof(float), 1, fp_dcd);
		coord[i].z = *dummyFloat;
	}

	fread(dummyInt, sizeof(int), 1, fp_dcd); // vectorByte
	
	/*
	printf("box_x :%lf, box_y : %lf, box_z : %lf\n",
			box[0], box[2], box[5]);
	*/

	/*
	for (int i = 0; i < 10; i++){
		printf("x : %f, y = %f, z = %f\n",
				coord[i].x, coord[i].y, coord[i].z);
	}
	*/

	return;
}

void initReader(){
	printf("Reading lammps_data ...\n");
	char *line = (char *)malloc(sizeof(char) * BUFFERSIZE);
	int id, type, ix, iy, iz;
	double x, y, z;

	for (int i = 0; i < 3; i++) fgets(line, BUFFERSIZE, fp_init); // IDK :D

	fscanf(fp_init, "%d %s", &numType, line);
	masses = (double *)malloc(sizeof(double) * numType);

	for (int i = 0; i < 7; i++) fgets(line, BUFFERSIZE, fp_init); // IDK :D

	for (int i = 0; i < numType; i++){
		fscanf(fp_init, "%d %lf", &type, &x);
		masses[type] = x;
	}

	for (int i = 0; i < 3; i++) fgets(line, BUFFERSIZE, fp_init); // IDK :D

	init = (particleInfo *)malloc(sizeof(int) * numAtom);

	for (int i = 0; i < numAtom; i++){
		fscanf(fp_init, "%d %d %lf %lf %lf %d %d %d",
				&id, &type, &x, &y, &z, &ix, &iy, &iz);
		// printf("%lf %lf %lf\n", x, y, z);
		init[i].id = id;
		init[i].type = type;
	}

	qsort(init, numAtom, sizeof(*init), compare); // sorting
	printf("lammps_data Reading Complete!\n");

	return;
}

int compare (const void * a, const void * b){
    return ( (*(particleInfo *)a).id - (*(particleInfo *)b).id );
}

void dcdHeader(){
	printf("Reading DCD Header ...\n");

	int ntitle;

	char *dummyChar = (char *)malloc(sizeof(char) * BUFFERSIZE);
	int *dummyInt = (int *)malloc(sizeof(int) * BUFFERSIZE);

	// Reading the header of dcd
	fseek(fp_dcd, 4, SEEK_CUR); // 84
	fseek(fp_dcd, 4, SEEK_CUR); // "CORD"

	fread(dummyInt , sizeof(int), 1, fp_dcd); // NFILE : # of snapshots in file
	numTraj = *dummyInt;

	fseek(fp_dcd, 96 - 12, SEEK_CUR); // IDK :D
	fread(dummyInt , sizeof(int), 1, fp_dcd); // ntitle
	ntitle = *dummyInt;

	fread(dummyChar, sizeof(char), ntitle * 80, fp_dcd); // title
	title = dummyChar;

	fseek(fp_dcd, 8, SEEK_CUR); // IDK :D
	fread(dummyInt , sizeof(int), 1, fp_dcd); // numAtom
	numAtom = *dummyInt;

	fseek(fp_dcd, 4, SEEK_CUR); // IDK :D

	// 112 + $ntitle * 80 bytes so far, theoretically.
	
	printf("\tnumTraj : %d\n", numTraj);
	printf("\tnumAtom : %d\n", numAtom);
	printf("DCD Header Reading complete!\n");

	return;
}
