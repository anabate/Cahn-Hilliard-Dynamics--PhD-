#include "cuda.h"
#include <stdio.h>
#include <math.h>
#include <curand.h>
#include <curand_kernel.h>
#include "gasdev.h"

#define HandleErrorWrapper( err ) (HandleError( err, __FILE__, __LINE__ ))


#define HandleNull( a ) {if (a == NULL) { \
                            printf( "Host memory failed in %s at line %d\n", \
                                    __FILE__, __LINE__ ); \
							waitKey();\
                            exit( EXIT_FAILURE );}}

// globals needed by the update routine
struct DataBlock {
    unsigned char   *output_bitmap;
    float           *dev_inSrc;
    float           *dev_outSrc;
	float           *dev_Energy;
	float			*TherF3D;
	float			*Pot2;
	float			*Pot3;

	curandState*    devStates;
    cudaEvent_t     start, stop;
    float           totalTime;
    float           frames;
};


// these exist on the GPU side
texture<float>  texIn;
texture<float>  texPot2;
texture<float>  texPot3;
texture<float>  texTherF3D;
texture<float>  texOut;
texture<float>  texE;

int saveDAT, savePLT, SaveEnergy, Nthread, hNz_count;

long Nsteps, NstepsClean,  
	 Vini, Metodo, errflg, 
	 hNx, hNy, hz_AB, z_AT, 
	 hNz, hNF, hNpot, Ntotal, 
	 bytes, idum, hHold; // or int vars

float *Fini, *X, *Y, *Z, *cFinif, *hPot2, *hPot3, *Energia;

float   htau, htauP, htauC, 
		ha, hf, hfP, hv, hu, hD, hDP, 
		hDx, hDy, hB, hBP, etaI, etaR, 
		hDT, hWd, hLa, hLb, hAmpl;
		
char ext[20];
errno_t io_error;

long *hNzs;

//initialization of the parameters
const char *FILE_PARAMS_BASENAME = "Input3D.dat";
//Save RandomMatrix
char FILE_RAND_BASENAME[FILENAME_MAX];
//Save PLT File Name
char FILE_PLT_BASENAME[FILENAME_MAX];
//Save DAT File Name
char FILE_DAT_BASENAME[FILENAME_MAX];
//Open a PLT file
char FILE_OPEN_BASENAME[FILENAME_MAX];
//Open an other PLT file
char FILE_OPEN_BASENAME_2[FILENAME_MAX];
// energy file 
char FILE_E_BASENAME[FILENAME_MAX];
//Variables to open a plt
char filename [ FILENAME_MAX ];
//Save PHI total DAT File Name
char FILE_PHI_TOTAL_DAT_BASENAME[FILENAME_MAX];
//Save PHI LAYERS DAT File Name
char FILE_PHI_DAT_BASENAME_CAPAS[FILENAME_MAX];

char cero [FILENAME_MAX];
char header1 [FILENAME_MAX];
char header2 [FILENAME_MAX];
char header3 [FILENAME_MAX];

//GPU constant memory
//Matrix Size
__constant__ int Nx;
__constant__ int Ny;
__constant__ int Nz;
__constant__ int z_AB;

//On off de congelar dinamica
__constant__ int Hold;

//1 --> agrega potencial 0 lo saca
__constant__ int Npot;

// Equivalent to Temperature (C: Clean)
__constant__ float tau;
__constant__ float tauC;
__constant__ float tauP;

//Map function parameters	
__constant__ float a;
__constant__ float v;
__constant__ float u;

//compositional assymetry
__constant__ float f;
__constant__ float fP;

//Diffusion Coefficient
__constant__ float D;
__constant__ float DP;

//Diffusion Coefficient in x
__constant__ float Dx;

//Diffusion Coefficient in y
__constant__ float Dy;

//Long range interaction
__constant__ float B;
__constant__ float BP;

//Noise level
__constant__ float eta;

//Time step size
__constant__ float DT;

//Potential parameters
__constant__ float Wd;
__constant__ float La;
__constant__ float Lb;
__constant__ float Ampl;


//************************************************************************************************************
//									FUNCTIONS CUDA
//*************************************************************************************************************
__global__ void setup_kernel ( curandState * state, unsigned long seed );
__device__ float generate( curandState* globalState, int ind );
__device__ float aver3D(texture<float> texMat, int x, int y, int z);
__device__ float aver3DPot(texture<float> texMat, texture<float> texC, int x, int y, int z);
__device__ float Derivx(texture<float> texMat, int x, int y, int z) ;
__device__ float Derivy(texture<float> texMat, int x, int y, int z) ;
__device__ float grad3D(texture<float> texMat, int x, int y, int z);
__global__ void gl_kernel1(	float *TherF3D, bool ruidoSwitch, bool dstOut);
__global__ void gl_kernel2( float *Fout, curandState* globalState, bool ruidoSwitch, bool dstOut);
__global__ void gl_kernel4( float *Energy, bool ruidoSwitch, bool dstOut);
void make_geo();
long calc_pot2();//float escala);
long calc_pot3();

//************************************************************************************************************
//									RUNNING and MEMORY_MANAGER
//*************************************************************************************************************
long init_mem();
void run_gl( DataBlock *d);
void anim_exit( DataBlock *d );
static void HandleError( cudaError_t err, const char *file, int line );

//************************************************************************************************************
//									FILE_MANAGER
//*************************************************************************************************************
errno_t get_plt(const char *fname);
errno_t get_plt_Alineado_BOTTOM(int hz_AB, const char *fname);
errno_t get_plt_Alineado_TOP(int hz_AB, const char *fname);
errno_t get_dat(char *fname);
errno_t get_param(const char *fname);
errno_t gen_rand(const char* fname, float *Fini);
errno_t gen_rand_Alineado_BOTTOM(int hz_AB, const char* fname, float *Fini);
errno_t save_plt(const char *basename, float *cFini, int q);
errno_t save_dat(const char *basename, float *cFini, int q);
errno_t save_dat_Energy(const char *basename, float EnergiaTotal, float FinalPhi, int q);
errno_t save_dat_Energy2(const char *basename, float EnergiaTotal, float FinalPhi, int k, int q);
errno_t save_dat_capas(const char *basename, float *energy, float *phi, int q);
void print_error(const char *msg, int err);
void waitKey();

//************************************************************************************************************
//									       MAIN
//*************************************************************************************************************

int main( void ) 
{
	printf("Initialization begin\n");

	//Loads the simulation parameters from file
	printf("Reading param file\n");
	if((io_error = get_param(FILE_PARAMS_BASENAME)) != 0)
		print_error("Error: opening param file", io_error);
	else
	{
		printf("OK \n");
	}
	
	printf("Verifying Nx is multiple of 16\n");
	if (!(hNx%16 == 0)) 
	{
		printf("Nx = %d is not a multiple of 16.\n", hNx); 
		getchar(); 
		exit(-1);
	}

	printf("Verifying Ny is multiple of 16\n");
	if (!(hNy%16 == 0)) 
	{
		printf("Ny = %d is not a multiple of 16\n", hNy); 
		getchar(); 
		exit(-1);
	}
	
	for(int i = 0; i < hNz_count; i++)
	{
		hNz = hNzs[i];
		
		// Calculate Ntotal
		printf("Processing Ntotal \n");
		Ntotal = hNx * hNy * hNz; 
		bytes = Ntotal*sizeof(float);

		// allocate memory for the needed arrays
		init_mem(); 
		printf("Verifying if the Potential is included..\n");
		if (hNpot == 1) 
		{
			calc_pot2();//(3.0f);
			calc_pot3();
			printf("The Potential is included..\n");
		} 
		else
		{
			printf("The Potential is NOT included..\n");
		}
		

	// if Vini = 0 I generate a random file, else I open other
	switch(Metodo)
	{
		case 0:
			if (Vini == 0)
			{
				printf("Generating random matrix...");
				if((io_error = gen_rand(FILE_RAND_BASENAME, Fini)) != 0)
				{
					print_error("Error: generating random matrix \n", io_error);
				}
				else
				{
					printf(" OK\n");
				}
			}
			break;
	
		case 1:	
			if(Vini != 0)
			{
				//// read Out file
				printf("Loading Initial file \n");
				if(strcmp(ext,"plt") == 0) 
				{
					if (Vini < 10) 
					{
						sprintf(cero, "0");
					}
					else 
					{
						sprintf(cero, "");
					}
					
					sprintf(filename, "%s%s%d.plt", FILE_OPEN_BASENAME, cero, Vini);
					printf("Loading %s... ",filename);
					errflg  =  get_plt(filename);
					if(errflg) 
					{
						print_error("Error: opening Initial .plt file", errflg);
					} 
					else 
					{
						printf("OK \n");
					}
				}

				if(strcmp(ext,"dat")==0) 
				{
					if (Vini<10) 
					{
						sprintf(cero, "0");
					} 
					else 
					{
						sprintf(cero, "");
					}
					
					sprintf(filename, "%s%s%d.dat", FILE_OPEN_BASENAME, cero, Vini);
					printf("Loading %s... ",filename);
					errflg  =  get_dat(filename);
					
					if(errflg) 
					{
						print_error("Error: opening Initial .dat file", errflg);
					} 
					else 
					{
						printf("OK \n");
					}
				}
			}
			break;
		
		case 2:
			if(Vini != 0)
				{ 
				printf("Loading Partial Initial file \n");
				if(strcmp(ext,"plt") == 0) 
				{
					if (Vini < 10) 
					{
						sprintf(cero, "0");
					}
					else 
					{
						sprintf(cero, "");
					}
					
					sprintf(filename, "%s%s%d.plt", FILE_OPEN_BASENAME, cero, Vini);
					printf("Loading %s... ",filename);
					errflg  =  get_plt_Alineado_BOTTOM(hz_AB, filename);
					
					if(errflg) 
					{
						print_error("Error: opening Initial .plt file", errflg);
					} 
					else 
					{
						printf("OK \n");
					}
				}
			

				printf("Generating random matrix...");
				if((io_error = gen_rand_Alineado_BOTTOM(hz_AB, FILE_RAND_BASENAME, Fini)) != 0)
				{
					print_error("Error: generating random matrix \n", io_error);
				}
				else
				{
					printf(" OK\n");
				}
			}
		break;	
		
		case 3:
			if(Vini != 0)
				{ 
				printf("Loading Partial Initial file \n");
				if(strcmp(ext,"plt") == 0) 
				{
					if (Vini < 10) 
					{
						sprintf(cero, "0");
					}
					else 
					{
						sprintf(cero, "");
					}
					
					sprintf(filename, "%s%s%d.plt", FILE_OPEN_BASENAME, cero, Vini);
					printf("Loading %s... ",filename);
					errflg  =  get_plt_Alineado_BOTTOM(hz_AB, filename);
			
					if(errflg) 
					{
						print_error("Error: opening Initial .plt file", errflg);
					} 
					else 
					{
						printf("OK \n");
					}
					
					sprintf(filename, "%s%s%d.plt", FILE_OPEN_BASENAME_2, cero, Vini);
					printf("Loading %s... ",filename);
					errflg  =  get_plt_Alineado_TOP(hz_AB, filename);
					
					if(errflg) 
					{
						print_error("Error: opening Initial .plt file", errflg);
					} 
					else 
					{
						printf("OK \n");
					}
				}
			}
		break;	
		
	}

		DataBlock   data;
		data.totalTime = 0;
		data.frames = 0;
		
		// random part
		if (!(etaR==0))
		{ 
			cudaMalloc( &data.devStates, Ntotal*sizeof( curandState ) );
		}
		// end of random
		
		HandleErrorWrapper( cudaEventCreate( &data.start ) );
		HandleErrorWrapper( cudaEventCreate( &data.stop ) );

		HandleErrorWrapper( cudaMalloc( (void**)&data.dev_inSrc, bytes ) );
		HandleErrorWrapper( cudaMalloc( (void**)&data.dev_outSrc, bytes ) );
		HandleErrorWrapper( cudaMalloc( (void**)&data.dev_Energy, bytes ) );
		HandleErrorWrapper( cudaMalloc( (void**)&data.TherF3D, bytes ) );
		HandleErrorWrapper( cudaMalloc( (void**)&data.Pot2, bytes ) ); 
		HandleErrorWrapper( cudaMalloc( (void**)&data.Pot3, bytes ) ); 
		
		HandleErrorWrapper( cudaBindTexture( NULL, texIn, data.dev_inSrc, bytes ) );
		HandleErrorWrapper( cudaBindTexture( NULL, texOut, data.dev_outSrc, bytes ) );
		HandleErrorWrapper( cudaBindTexture( NULL, texE, data.dev_Energy, bytes ) );
		HandleErrorWrapper( cudaBindTexture( NULL, texTherF3D, data.TherF3D, bytes ) );
		HandleErrorWrapper( cudaBindTexture( NULL, texPot2, data.Pot2, bytes ) );
		HandleErrorWrapper( cudaBindTexture( NULL, texPot3, data.Pot3, bytes ) );
			
		// CUDA memory copy
		HandleErrorWrapper( cudaMemcpy( data.dev_inSrc, Fini, bytes, cudaMemcpyHostToDevice ) );
		HandleErrorWrapper( cudaMemcpy( data.TherF3D, Fini, bytes, cudaMemcpyHostToDevice ) );
		HandleErrorWrapper( cudaMemcpy( data.Pot2, hPot2, bytes, cudaMemcpyHostToDevice ) );
		HandleErrorWrapper( cudaMemcpy( data.Pot3, hPot3, bytes, cudaMemcpyHostToDevice ) );
		
		cudaMemcpyToSymbol(Nx, &hNx, sizeof(int) );
		cudaMemcpyToSymbol(Ny, &hNy, sizeof(int) );
		cudaMemcpyToSymbol(Nz, &hNz, sizeof(int) );
		cudaMemcpyToSymbol(z_AB, &hz_AB, sizeof(int) );
		
		cudaMemcpyToSymbol(Npot, &hNpot, sizeof(int) );
		cudaMemcpyToSymbol(Hold, &hHold, sizeof(int) );
		cudaMemcpyToSymbol(tau, &htau, sizeof(float) );
		cudaMemcpyToSymbol(tauP, &htauP, sizeof(float) );
		cudaMemcpyToSymbol(tauC, &htauC, sizeof(float) );
		cudaMemcpyToSymbol(a, &ha, sizeof(float) );
		cudaMemcpyToSymbol(f, &hf, sizeof(float) );
		cudaMemcpyToSymbol(fP, &hfP, sizeof(float) );
		cudaMemcpyToSymbol(v, &hv, sizeof(float) );
		cudaMemcpyToSymbol(u, &hu, sizeof(float) );
		cudaMemcpyToSymbol(D, &hD, sizeof(float) );
		cudaMemcpyToSymbol(DP, &hDP, sizeof(float) );
		cudaMemcpyToSymbol(Dx, &hDx, sizeof(float) );
		cudaMemcpyToSymbol(Dy, &hDy, sizeof(float) );
		cudaMemcpyToSymbol(B, &hB, sizeof(float) );
		cudaMemcpyToSymbol(BP, &hBP, sizeof(float) );
		cudaMemcpyToSymbol(eta, &etaR, sizeof(float) );
		cudaMemcpyToSymbol(DT, &hDT, sizeof(float) );
		cudaMemcpyToSymbol(Wd, &hWd, sizeof(float) );
		cudaMemcpyToSymbol(La, &hLa, sizeof(float) );
		cudaMemcpyToSymbol(Lb, &hLb, sizeof(float) );
		cudaMemcpyToSymbol(Ampl, &hAmpl, sizeof(float) );


		run_gl(&data);
		anim_exit( &data );
	}
	
	waitKey();
}

//************************************************************************************************************
//									FUNCTIONS CUDA
//*************************************************************************************************************

__global__ void setup_kernel ( curandState * state, unsigned long seed )
{
	int x = threadIdx.x + blockIdx.x * blockDim.x; 
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    int irand = x + y * blockDim.x * gridDim.x;
	curand_init ( seed, irand, 0, &state[irand] );
}

__device__ float generate( curandState* globalState, int ind ) 
{
    curandState localState = globalState[ind];
	float RANDOM = curand_uniform( &localState );
	float RANDOM2 = eta*(1.0f-2.0f*RANDOM);
    globalState[ind] = localState;
    return RANDOM2;
}

__device__ float aver3D(texture<float> texMat, int x, int y, int z) 
{
	float average;
	
	int ixl = x - 1;
	int ixr = x + 1;
	
	if (x == 0)
	{
		ixl=(Nx-1);
	}
	
	if (x == Nx-1) 
	{
		ixr=0; 
	}
	int iyl = y - 1;
	int iyr = y + 1;
	if (y == 0)
	{
		iyl=(Ny-1);
	}
	if (y == Ny-1) 
	{
		iyr=0; 
	}
	int izl = z - 1;
	int izr = z + 1;
	
	if (z == 0)
	{
		izl=(Nz-1);
	}
	
	if (z == Nz-1)
	{
		izr=0; 
	}
		
	float Xl, Xr, Yl, Yr, Zl, Zr, Xlt, Xrt, Xrb, Xlb, Zrl, Zrr, Zrt, Zrb, Zll, Zlr, Zlt, Zlb, Zrlt, Zrrt, Zrrb, Zrlb, Zllt, Zlrt, Zlrb, Zllb;
	
	int Nyy, Nyr, Nyl, Nzz, Nzr, Nzl;
	Nyy = y * Nx; 
	Nyr = iyr * Nx; 
	Nyl = iyl * Nx;
	Nzz = z * Nx * Ny; 
	Nzr = izr * Nx * Ny; 
	Nzl = izl * Nx * Ny;
			
	// first neighbors (6)
	Xl = tex1Dfetch(texMat,ixl + Nyy + Nzz);
	Xr = tex1Dfetch(texMat,ixr + Nyy + Nzz);
	Yl = tex1Dfetch(texMat,x + Nyl + Nzz);
	Yr = tex1Dfetch(texMat,x + Nyr + Nzz);
	Zl = tex1Dfetch(texMat,x + Nyy + Nzl);
	Zr = tex1Dfetch(texMat,x + Nyy + Nzr);
	
	// second neighbors (12)
	Xlt = tex1Dfetch(texMat,ixl + Nyr + Nzz);
	Xrt = tex1Dfetch(texMat,ixr + Nyr + Nzz);
	Xrb = tex1Dfetch(texMat,ixr + Nyl + Nzz);
	Xlb = tex1Dfetch(texMat,ixl + Nyl + Nzz);
	
	Zrl = tex1Dfetch(texMat,ixl + Nyy + Nzr);
	Zrr = tex1Dfetch(texMat,ixr + Nyy + Nzr);
	Zrt = tex1Dfetch(texMat,x + Nyr + Nzr);
	Zrb = tex1Dfetch(texMat,x + Nyl + Nzr);
	 
	Zll = tex1Dfetch(texMat,ixl + Nyy + Nzl);
	Zlr = tex1Dfetch(texMat,ixr + Nyy + Nzl);
	Zlt = tex1Dfetch(texMat,x + Nyr + Nzl);
	Zlb = tex1Dfetch(texMat,x + Nyl + Nzl);
	
	// third neighbors (8)
	Zrlt = tex1Dfetch(texMat,ixl + Nyr + Nzr);
	Zrrt = tex1Dfetch(texMat,ixr + Nyr + Nzr);
	Zrrb = tex1Dfetch(texMat,ixr + Nyl + Nzr);
	Zrlb = tex1Dfetch(texMat,ixl + Nyl + Nzr);
	
	Zllt = tex1Dfetch(texMat,ixl + Nyr + Nzl);
	Zlrt = tex1Dfetch(texMat,ixr + Nyr + Nzl);
	Zlrb = tex1Dfetch(texMat,ixr + Nyl + Nzl);
	Zllb = tex1Dfetch(texMat,ixl + Nyl + Nzl);
	
	float c_1 = 6.0f/80.0f;
	float c_2 = 3.0f/80.0f;
	float c_3 = 1.0f/80.0f;
	
	average = c_1*(Xl + Xr + Yl + Yr + Zl + Zr) + 
			  c_2*(Xlt + Xrt + Xrb + Xlb + Zrl + Zrr + Zrt + Zrb + Zll + Zlr + Zlt + Zlb) + 
			  c_3*(Zrlt + Zrrt + Zrrb + Zrlb + Zllt + Zlrt + Zlrb + Zllb);
	
	return average;
}


__device__ float aver3DPot(texture<float> texMat, texture<float> texC, int x, int y, int z) 
{
	float average;
	
	int ixl = x - 1;
	int ixr = x + 1;
	
	if (x == 0)
	{
		ixl=(Nx-1);
	}
	
	if (x == Nx-1) 
	{
		ixr=0; 
	}
	int iyl = y - 1;
	int iyr = y + 1;
	if (y == 0)    iyl=(Ny-1);
	if (y == Ny-1) iyr=0; 
	
	int izl = z - 1;
	int izr = z + 1;
	
	if (z == 0)
	{
		izl=(Nz-1);
	}
	
	if (z == Nz-1)
	{
		izr=0; 
	}
		
	float Xl, Xr, Yl, Yr, Zl, Zr, Xlt, Xrt, Xrb, Xlb, Zrl, Zrr, Zrt, Zrb, Zll, Zlr, Zlt, Zlb, Zrlt, Zrrt, Zrrb, Zrlb, Zllt, Zlrt, Zlrb, Zllb;
	
	int Nyy, Nyr, Nyl, Nzz, Nzr, Nzl;
	Nyy = y * Nx; 
	Nyr = iyr * Nx; 
	Nyl = iyl * Nx;
	Nzz = z * Nx * Ny; 
	Nzr = izr * Nx * Ny; 
	Nzl = izl * Nx * Ny;
			
	// first neighbors (6)
	Xl = tex1Dfetch(texMat,ixl + Nyy + Nzz) *  tex1Dfetch(texC,ixl + Nyy + Nzz);
	Xr = tex1Dfetch(texMat,ixr + Nyy + Nzz) *  tex1Dfetch(texC,ixr + Nyy + Nzz);
	Yl = tex1Dfetch(texMat,x   + Nyl + Nzz) *  tex1Dfetch(texC,x   + Nyl + Nzz);
	Yr = tex1Dfetch(texMat,x   + Nyr + Nzz) *  tex1Dfetch(texC,x   + Nyr + Nzz);
	Zl = tex1Dfetch(texMat,x   + Nyy + Nzl) *  tex1Dfetch(texC,x   + Nyy + Nzl);
	Zr = tex1Dfetch(texMat,x   + Nyy + Nzr) *  tex1Dfetch(texC,x   + Nyy + Nzr);
	
	// secneighbors (12)
	Xlt = tex1Dfetch(texMat,ixl + Nyr + Nzz) * 	tex1Dfetch(texC,ixl + Nyr + Nzz);
	Xrt = tex1Dfetch(texMat,ixr + Nyr + Nzz) *  tex1Dfetch(texC,ixr + Nyr + Nzz);
	Xrb = tex1Dfetch(texMat,ixr + Nyl + Nzz) *  tex1Dfetch(texC,ixr + Nyl + Nzz);
	Xlb = tex1Dfetch(texMat,ixl + Nyl + Nzz) *  tex1Dfetch(texC,ixl + Nyl + Nzz);
	
	Zrl = tex1Dfetch(texMat,ixl + Nyy + Nzr) *  tex1Dfetch(texC,ixl + Nyy + Nzr);
	Zrr = tex1Dfetch(texMat,ixr + Nyy + Nzr) *  tex1Dfetch(texC,ixr + Nyy + Nzr);
	Zrt = tex1Dfetch(texMat,x   + Nyr + Nzr) *  tex1Dfetch(texC,x   + Nyr + Nzr);;
	Zrb = tex1Dfetch(texMat,x   + Nyl + Nzr) *  tex1Dfetch(texC,x   + Nyl + Nzr);
	 	  
	Zll = tex1Dfetch(texMat,ixl + Nyy + Nzl) *  tex1Dfetch(texC,ixl + Nyy + Nzl);
	Zlr = tex1Dfetch(texMat,ixr + Nyy + Nzl) *  tex1Dfetch(texC,ixr + Nyy + Nzl);
	Zlt = tex1Dfetch(texMat,x   + Nyr + Nzl) *  tex1Dfetch(texC,x   + Nyr + Nzl);;
	Zlb = tex1Dfetch(texMat,x   + Nyl + Nzl) *  tex1Dfetch(texC,x   + Nyl + Nzl);
	
		 
	// thieighbors (8)
	Zrlt = tex1Dfetch(texMat,ixl + Nyr + Nzr) * tex1Dfetch(texC,ixl + Nyr + Nzr);
	Zrrt = tex1Dfetch(texMat,ixr + Nyr + Nzr) * tex1Dfetch(texC,ixr + Nyr + Nzr);
	Zrrb = tex1Dfetch(texMat,ixr + Nyl + Nzr) * tex1Dfetch(texC,ixr + Nyl + Nzr);
	Zrlb = tex1Dfetch(texMat,ixl + Nyl + Nzr) * tex1Dfetch(texC,ixl + Nyl + Nzr);
		  													  
	Zllt = tex1Dfetch(texMat,ixl + Nyr + Nzl) * tex1Dfetch(texC,ixl + Nyr + Nzl);
	Zlrt = tex1Dfetch(texMat,ixr + Nyr + Nzl) * tex1Dfetch(texC,ixr + Nyr + Nzl);
	Zlrb = tex1Dfetch(texMat,ixr + Nyl + Nzl) * tex1Dfetch(texC,ixr + Nyl + Nzl);
	Zllb = tex1Dfetch(texMat,ixl + Nyl + Nzl) * tex1Dfetch(texC,ixl + Nyl + Nzl);

	
	float c_1 = 6.0f/80.0f;
	float c_2 = 3.0f/80.0f;
	float c_3 = 1.0f/80.0f;
	
	average =    c_1*(Xl + Xr + Yl + Yr + Zl + Zr) + 
		    	 c_2*(Xlt + Xrt + Xrb + Xlb + Zrl + Zrr + Zrt + Zrb + Zll + Zlr + Zlt + Zlb) + 
				 c_3*(Zrlt + Zrrt + Zrrb + Zrlb + Zllt + Zlrt + Zlrb + Zllb);
	
	return average;
	
}


__device__ float grad3D(texture<float> texMat, int x, int y, int z) 
{
	float gradient, gradient2, dx, dy, dz, c, delta;
	
	int ixl = x - 1;
	int ixr = x + 1;
	
	if (x == 0)
	{
		ixl=(Nx-1);
	}
	
	if (x == Nx-1) 
	{
		ixr=0; 
	}
	int iyl = y - 1;
	int iyr = y + 1;
	if (y == 0)
	{
		iyl=(Ny-1);
	}
	if (y == Ny-1) 
	{
		iyr=0; 
	}
	int izl = z - 1;
	int izr = z + 1;
	
	if (z == 0)
	{
		izl=(Nz-1);
	}
	
	if (z == Nz-1)
	{
		izr=0; 
	}


	float Xr, Yr, Zr, Xl, Yl, Zl, Xlt;
	
	int Nyy, Nyr, Nyl, Nzz, Nzr, Nzl;
	
	delta = powf((40.0f/11.0f),(0.5f));

	Nyy = y * Nx; 
	Nyr = iyr * Nx; 
	Nyl = iyl * Nx;
	Nzz = z * Nx * Ny; 
	Nzr = izr * Nx * Ny; 
	Nzl = izl * Nx * Ny;
	
	c = tex1Dfetch(texMat, x + Nyy + Nzz);
	
	
	Xl = tex1Dfetch(texMat,ixl + Nyy + Nzz);
	Xr = tex1Dfetch(texMat,ixr + Nyy + Nzz);
	
	Yl = tex1Dfetch(texMat,x + Nyl + Nzz);
	Yr = tex1Dfetch(texMat,x + Nyr + Nzz);
	
	Zl = tex1Dfetch(texMat,x + Nyy + Nzl);
	Zr = tex1Dfetch(texMat,x + Nyy + Nzr);

	
	dx = (Xr - Xl)/(2.0f*delta);
	dy = (Yr - Yl)/(2.0f*delta);
	dz = (Zr - Zl)/(2.0f*delta);
	
	gradient2 = powf(dx,2) + powf(dy,2) + powf(dz,2);
	
	gradient =  powf(gradient2,0.5f);
	
	return gradient;
	
}

__device__ float Derivx(texture<float> texMat, int x, int y, int z) 
{
	float d2x;
	
	int ixl = x - 1;
	int ixr = x + 1;
	
	if (x == 0)
	{
		ixl=(Nx-1);
	}
		if (x == Nx-1) 
	{
		ixr=0; 
	}
			
	float Xl, Xr;
	float c_1 = 1.0f/2.0f;
	int Nyy, Nzz;
	
	Nyy = y * Nx; 
	Nzz = z * Nx * Ny; 
					
	Xl = tex1Dfetch(texMat,ixl + Nyy + Nzz);
	Xr = tex1Dfetch(texMat,ixr + Nyy + Nzz);
		
	d2x = c_1*(Xr + Xl) ;
	
	return d2x;
}

__device__ float Derivy(texture<float> texMat, int x, int y, int z) 
{
	float d2y;
	
	int iyl = y - 1;
	int iyr = y + 1;
	
	if (y == 0)
	{
		iyl=(Ny-1);
	}
	if (y == Ny-1)
	{
		iyr=0; 
	}
	
		
	float Yl, Yr;
	
	int Nyr, Nyl, Nzz;
	
	Nyr = iyr * Nx; 
	Nyl = iyl * Nx;
	Nzz = z * Nx * Ny; 
	
			
	Yl = tex1Dfetch(texMat,x + Nyl + Nzz);
	Yr = tex1Dfetch(texMat,x + Nyr + Nzz);
	
	float c_1 = 1.0f/2.0f;
	
	d2y = c_1*(Yr + Yl);
	
	return d2y;
	
}

__global__ void gl_kernel1(	float *TherF3D, bool ruidoSwitch, bool dstOut)
{
    // map from threadIdx/BlockIdx to pixel position
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
	int z = threadIdx.z + blockIdx.z * blockDim.z;
	
	int i = x + y * Nx + z * Nx * Ny;
	
	float c, average, tauL, d2x, d2y;
			
	if (x < Nx && y < Ny && z < Nz) 
	{
		//ENTRO SI: hold es 0 entro//  hold es 1 y z mayor a z_AB 
		
		switch(Hold)
		{
			case 0: //MONOCAPA
				
				float c_1 = 1.0f - 2.0f*f;
				if (dstOut) 
				{
					c = tex1Dfetch(texIn, i); // center value of Fini
					average = aver3D(texIn, x, y, z);
					d2x = Derivx(texIn, x, y, z);
					d2y = Derivy(texIn, x, y, z);
				} 
				else 
				{// center value of Fini
					c = tex1Dfetch(texOut, i); 
					average = aver3D(texOut, x, y, z);
					d2x = Derivx(texOut, x, y, z);
					d2y = Derivy(texOut, x, y, z);
				}
			
			  
				//// Solving Ginzburg-Landau Equation
				//TherF3D // TherF3D=Fm-D*(Aver3D(Fini)-Fini) // Aver3D := laplacian // Fm=MapF3D(Fini)
											
				TherF3D[i] = -(tau - a * powf(c_1,2)) * c
									+ v * c_1 *powf(c,2)
									+ u * powf(c,3)  
									- D * (average - c) 
									- Dx * (d2x - c) 
									- Dy * (d2y - c) ;



			case 1: // BIACPA: PARTE DE ABAJO ES FIJA

				if (z > z_AB)
				{
					float c_1 = 1.0f - 2.0f*fP;

					if (dstOut) 
					{
						c = tex1Dfetch(texIn, i); // center value of Fini
						average = aver3D(texIn, x, y, z);
						d2x = Derivx(texIn, x, y, z);
						d2y = Derivy(texIn, x, y, z);
					} 
					else 
					{
						c = tex1Dfetch(texOut, i); // center value of Fini
						average = aver3D(texOut, x, y, z);
						d2x = Derivx(texOut, x, y, z);
						d2y = Derivy(texOut, x, y, z);
					}
			
			//// Solving Ginzburg-Landau Equation
			//TherF3D // TherF3D=Fm-D*(Aver3D(Fini)-Fini) // Aver3D := laplacian // Fm=MapF3D(Fini)
			
								
					TherF3D[i]=-(tauP - a * powf(c_1,2)) * c
										  +	v * c_1 *powf(c,2)
										  +	u * powf(c,3)  
										  -	(DP) * (average - c) 
										  - Dx * (d2x - c) 
										  - Dy * (d2y - c) ;							
				}

				else
				{
					float c_1 = 1.0f - 2.0f*f;
					if (dstOut) 
					{
						c = tex1Dfetch(texIn, i); // center value of Fini
					} 
					else 
					{// center value of Fini
						c = tex1Dfetch(texOut, i); 						
					}
				
						TherF3D[i] = c;		//	 ABAJO FIJO ES SOLO IGUAL A C		
				}


				case 2: // BIACPA: PARTE DE ABAJO EVOLUCIONA TAMBIEN
					
				if (z > z_AB)
				{
					float c_1 = 1.0f - 2.0f*fP;

					if (dstOut) 
					{
						c = tex1Dfetch(texIn, i); // center value of Fini
						average = aver3D(texIn, x, y, z);
						d2x = Derivx(texIn, x, y, z);
						d2y = Derivy(texIn, x, y, z);
					} 
					else 
					{
						c = tex1Dfetch(texOut, i); // center value of Fini
						average = aver3D(texOut, x, y, z);
						d2x = Derivx(texOut, x, y, z);
						d2y = Derivy(texOut, x, y, z);
					}
			
			//// Solving Ginzburg-Landau Equation
			//TherF3D // TherF3D=Fm-D*(Aver3D(Fini)-Fini) // Aver3D := laplacian // Fm=MapF3D(Fini)
			
								
						TherF3D[i]=-(tauP - a * powf(c_1,2)) * c
										  +	v * c_1 *powf(c,2)
										  +	u * powf(c,3)  
										  -	(DP) * (average - c) 
										  - Dx * (d2x - c) 
										  - Dy * (d2y - c) ;
								
				}
				else
				{
					float c_1 = 1.0f - 2.0f*f;
					if (dstOut) 
					{
						c = tex1Dfetch(texIn, i); // center value of Fini
						average = aver3D(texIn, x, y, z);
						d2x = Derivx(texIn, x, y, z);
						d2y = Derivy(texIn, x, y, z);
					} 
					else 
					{// center value of Fini
						c = tex1Dfetch(texOut, i); 
						average = aver3D(texOut, x, y, z);
						d2x = Derivx(texOut, x, y, z);
						d2y = Derivy(texOut, x, y, z);
					}
				
			
			  
					//// Solving Ginzburg-Landau Equation
					//TherF3D // TherF3D=Fm-D*(Aver3D(Fini)-Fini) // Aver3D := laplacian // Fm=MapF3D(Fini)
											
						TherF3D[i] = -(tau - a * powf(c_1,2)) * c
										+ v * c_1 *powf(c,2)
										+ u * powf(c,3)  
										- D * (average - c) 
										- Dx * (d2x - c) 
										- Dy * (d2y - c) ;
				}
		}
	}
} 

__global__ void gl_kernel2( float *Fout, curandState* globalState, bool ruidoSwitch, bool dstOut) 
{
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
	int z = threadIdx.z + blockIdx.z * blockDim.z;
	
	int irand = x + y * blockDim.x * gridDim.x;
	
	int i = x + y * Nx + z * Nx * Ny;
	
	float c, TF, TFav, CenterPot, CenterPot2, AverPot, AverPot2;
		
	
	if (x < Nx && y < Ny && z < Nz) 
	{			
		switch(Hold)
		{
			case 0: //MONOCAPA

				if (dstOut) 
				{
					c = tex1Dfetch(texIn, i); 
					AverPot = aver3DPot(texPot2, texIn, x, y, z); 
				
				} 
				else 
				{
					c = tex1Dfetch(texOut, i);
					AverPot = aver3DPot(texPot2, texOut, x, y, z); 
				}	

				CenterPot = tex1Dfetch(texPot2, i);// center value of texPot2c

				TF = tex1Dfetch(texTherF3D, i); 

				TFav = aver3D(texTherF3D, x, y, z);
		   
				AverPot2 = aver3D(texPot3, x, y, z);
				CenterPot2 = tex1Dfetch(texPot3, i);// center value of texPot2c

				float noise = 0.0f;
				if(ruidoSwitch) 
				{
					if (!(eta==0.0f)) 
					{
						noise=generate(globalState, irand);
					}
				}
				if (Npot) 
				{ 
					Fout[i] = c + DT * (TFav - TF - B * c) +  DT * (AverPot - CenterPot * c ) + DT * (AverPot2 - CenterPot2) +  noise ;
    			}
				else
				{ 	
					Fout[i] = c + DT * (TFav - TF - (B) * c)  +   noise ;
				}


			case 1: // BIACPA: PARTE DE ABAJO ES FIJA

			if (z > z_AB)
			{
			
				if (dstOut) 
				{
					c = tex1Dfetch(texIn, i); 
					AverPot = aver3DPot(texPot2, texIn, x, y, z); 
				} 
				else 
				{
					c = tex1Dfetch(texOut, i);
					AverPot = aver3DPot(texPot2, texOut, x, y, z); 
				}	

				TF = tex1Dfetch(texTherF3D, i); 

				TFav = aver3D(texTherF3D, x, y, z);

				CenterPot = tex1Dfetch(texPot2, i); // center value of texPot2
			
				AverPot2 = aver3D(texPot3, x, y, z);
				CenterPot2 = tex1Dfetch(texPot3, i); // center value of texPot3
				
			
				float noise = 0.0f;
				if(ruidoSwitch) 
				{
					if (!(eta==0.0f)) 
					{
						noise=generate(globalState, irand);
					}
				}
				
				if (Npot) 
				{ 
					Fout[i] = c + DT * (TFav - TF - (BP) * c)  +  DT * (AverPot - CenterPot * c ) + DT * (AverPot2 - CenterPot2) + noise;
				}
				else
				{ 
					Fout[i] = c + DT * (TFav - TF - (BP) * c)  +  noise;
				}

			}

			else
			{
				if (dstOut) 
				{
					c = tex1Dfetch(texIn, i); 
					AverPot = aver3DPot(texPot2, texIn, x, y, z); 
				
				} 
				else 
				{
					c = tex1Dfetch(texOut, i);
					AverPot = aver3DPot(texPot2, texOut, x, y, z); 
				}	

				CenterPot = tex1Dfetch(texPot2, i);// center value of texPot2c  
				AverPot2 = aver3D(texPot3, x, y, z);
				CenterPot2 = tex1Dfetch(texPot3, i);// center value of texPot2c

				float noise = 0.0f;
				if(ruidoSwitch) 
				{
					if (!(eta==0.0f)) 
					{
						noise=generate(globalState, irand);
					}
				}
				if (Npot) 
				{ 
						 Fout[i] = c + DT * (AverPot - CenterPot * c ) + DT * (AverPot2 - CenterPot2);  //ABAJO FIJO ES ESTO
				}
				else
				{ 	
						Fout[i] = c ;// ABAJO FIJO ES ESTO
				}
		}


			case 2: // BIACPA: PARTE DE ABAJO TIENE DINAMICA

			if (z > z_AB)
			{
			
				if (dstOut) 
				{
					c = tex1Dfetch(texIn, i); 
					AverPot = aver3DPot(texPot2, texIn, x, y, z); 
				} 
				else 
				{
					c = tex1Dfetch(texOut, i);
					AverPot = aver3DPot(texPot2, texOut, x, y, z); 
				}	

				TF = tex1Dfetch(texTherF3D, i); 

				TFav = aver3D(texTherF3D, x, y, z);

				CenterPot = tex1Dfetch(texPot2, i); // center value of texPot2
				AverPot2 = aver3D(texPot3, x, y, z);
				CenterPot2 = tex1Dfetch(texPot3, i); // center value of texPot3
				
			
				float noise = 0.0f;
				if(ruidoSwitch) 
				{
					if (!(eta==0.0f)) 
					{
						noise=generate(globalState, irand);
					}
				}
				
				if (Npot) 
				{ 
					Fout[i] = c + DT * (TFav - TF - (BP) * c)  +  DT * (AverPot - CenterPot * c ) + DT * (AverPot2 - CenterPot2) + noise;
				}
				else
				{ 
					Fout[i] = c + DT * (TFav - TF - (BP) * c)  +  noise;
				}

			}

			else
			{
				if (dstOut) 
				{
					c = tex1Dfetch(texIn, i); 
					AverPot = aver3DPot(texPot2, texIn, x, y, z); 
				
				} 
				else 
				{
					c = tex1Dfetch(texOut, i);
					AverPot = aver3DPot(texPot2, texOut, x, y, z); 
				}	

				CenterPot = tex1Dfetch(texPot2, i);// center value of texPot2c

				TF = tex1Dfetch(texTherF3D, i); 

				TFav = aver3D(texTherF3D, x, y, z);
		   
				AverPot2 = aver3D(texPot3, x, y, z);
				CenterPot2 = tex1Dfetch(texPot3, i);// center value of texPot2c

				float noise = 0.0f;
				if(ruidoSwitch) 
				{
					if (!(eta==0.0f)) 
					{
						noise=generate(globalState, irand);
					}
				}
				if (Npot) 
				{ 
					Fout[i] = c + DT * (TFav - TF - B * c) +  DT * (AverPot - CenterPot * c ) + DT * (AverPot2 - CenterPot2) +  noise ;
    			}
				else
				{ 	
					Fout[i] = c + DT * (TFav - TF - (B) * c)  +   noise ;
				}
			}
		}
	}
}

			
__global__ void gl_kernel4( float *Energy, bool ruidoSwitch, bool dstOut) 
{
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
	int z = threadIdx.z + blockIdx.z * blockDim.z;
	
	int irand = x + y * blockDim.x * gridDim.x;
	
	int i = x + y * Nx + z * Nx * Ny;
	
	float tauL, CenterLaplacian3, AverLaplacian3, GradLaplacian3, Laplaciano;

				
	if (x<Nx && y<Ny && z<Nz) 
	{	

		switch(Hold)
		{
			case 0: //MONOCAPA
				
				float c_1 = 0.5f*(powf(D,3/2)/powf(B,0.5f));
				float c_2 = - D;
				float c_3 = (0.5f)*(-tau + a*powf((1.0f-2.0f*f),2) + 3 * powf(B*D, 0.5f)  );
				float c_4 = (v/3)*(1.0f-2.0f*fP);
				float c_5 = (u/4);

				if (dstOut) 
				{
					CenterLaplacian3 = tex1Dfetch(texIn, i); 
					AverLaplacian3 = aver3D(texIn, x, y, z);
					GradLaplacian3 = grad3D(texIn, x, y, z);
				} 
				else 
				{
					CenterLaplacian3 = tex1Dfetch(texOut, i); 
					AverLaplacian3 = aver3D(texOut, x, y, z);
					GradLaplacian3 = grad3D(texOut, x, y, z);
				}
				
				if (ruidoSwitch) 
				{
					tauL = tau; 
				} 
				else 
				{
					tauL = tauC;
				}
				
				Laplaciano = AverLaplacian3 - CenterLaplacian3;
			
				Energy[i] =     	c_1 * powf(Laplaciano,2) +
									c_2 * powf(GradLaplacian3,2) + 
									c_3 * powf(CenterLaplacian3,2) + 
									c_4 * powf(CenterLaplacian3,3) ;//+
									c_5 * powf(CenterLaplacian3,4);		
		
			case 1: // BIACPA: PARTE DE ABAJO ES FIJA

				if (z > z_AB)
				{
					float c_1 = 0.5f*(powf(DP,3/2)/powf(BP,0.5f));
					float c_2 = - DP;
					float c_3 = (0.5f)*(-tauP + a*powf((1.0f-2.0f*fP),2) + 3 * powf(BP*DP, 0.5f)  );
					float c_4 = (v/3)*(1.0f-2.0f*fP);
					float c_5 = (u/4);

					if (dstOut) 
					{
						CenterLaplacian3 = tex1Dfetch(texIn, i); 
						AverLaplacian3 = aver3D(texIn, x, y, z);
						GradLaplacian3 = grad3D(texIn, x, y, z);
					} 
					else 
					{
						CenterLaplacian3 = tex1Dfetch(texOut, i); 
						AverLaplacian3 = aver3D(texOut, x, y, z);
						GradLaplacian3 = grad3D(texOut, x, y, z);
					}
				
					if (ruidoSwitch) 
					{
						tauL = tau; 
					} 
					else 
					{
						tauL = tauC;
					}
				
						Laplaciano = AverLaplacian3 - CenterLaplacian3;
			
						Energy[i] =     	c_1 * powf(Laplaciano,2) +
											c_2 * powf(GradLaplacian3,2) + 
											c_3 * powf(CenterLaplacian3,2) + 
											c_4 * powf(CenterLaplacian3,3);// +
											c_5 * powf(CenterLaplacian3,4);		
				}
		
				else
				{

					float c_1 = 0.5f*(powf(D,3/2)/powf(B,0.5f));
					float c_2 = - D;
					float c_3 = (0.5f)*(-tau + a*powf((1.0f-2.0f*f),2) + 3 * powf(B*D, 0.5f)  );
					float c_4 = (v/3)*(1.0f-2.0f*fP);
					float c_5 = (u/4);

					if (dstOut) 
					{
						CenterLaplacian3 = tex1Dfetch(texIn, i); 
						AverLaplacian3 = aver3D(texIn, x, y, z);
						GradLaplacian3 = grad3D(texIn, x, y, z);
					} 
					else 
					{
						CenterLaplacian3 = tex1Dfetch(texOut, i); 
						AverLaplacian3 = aver3D(texOut, x, y, z);
						GradLaplacian3 = grad3D(texOut, x, y, z);
					}
				
					if (ruidoSwitch) 
					{
						tauL = tau; 
					} 
					else 
					{
						tauL = tauC;
					}
				
					Laplaciano = AverLaplacian3 - CenterLaplacian3;
			
					Energy[i] =     	c_1 * powf(Laplaciano,2) +
										c_2 * powf(GradLaplacian3,2) + 
										c_3 * powf(CenterLaplacian3,2) + 
										c_4 * powf(CenterLaplacian3,3) +
										c_5 * powf(CenterLaplacian3,4);		
				}

			case 2: // BIACPA: PARTE DE ABAJO TIENE DINAMICA

				if (z > z_AB)
				{
					float c_1 = 0.5f*(powf(DP,3/2)/powf(BP,0.5f));
					float c_2 = - DP;
					float c_3 = (0.5f)*(-tauP + a*powf((1.0f-2.0f*fP),2) + 3 * powf(BP*DP, 0.5f)  );
					float c_4 = (v/3)*(1.0f-2.0f*fP);
					float c_5 = (u/4);

					if (dstOut) 
					{
						CenterLaplacian3 = tex1Dfetch(texIn, i); 
						AverLaplacian3 = aver3D(texIn, x, y, z);
						GradLaplacian3 = grad3D(texIn, x, y, z);
					} 
					else 
					{
						CenterLaplacian3 = tex1Dfetch(texOut, i); 
						AverLaplacian3 = aver3D(texOut, x, y, z);
						GradLaplacian3 = grad3D(texOut, x, y, z);
					}
				
					if (ruidoSwitch) 
					{
						tauL = tau; 
					} 
					else 
					{
						tauL = tauC;
					}
				
						Laplaciano = AverLaplacian3 - CenterLaplacian3;
			
						Energy[i] =     	c_1 * powf(Laplaciano,2) +
											c_2 * powf(GradLaplacian3,2) + 
											c_3 * powf(CenterLaplacian3,2) + //este hay q sacar
											c_4 * powf(CenterLaplacian3,3);// +
											c_5 * powf(CenterLaplacian3,4);		
				}
		
				else
				{

					float c_1 = 0.5f*(powf(D,3/2)/powf(B,0.5f));
					float c_2 = - D;
					float c_3 = -(0.5f)*(tau - a*powf((1.0f-2.0f*f),2) + 3 * powf(B*D, 0.5f)  );
					float c_4 = (v/3)*(1.0f-2.0f*fP);
					float c_5 = (u/4);

					if (dstOut) 
					{
						CenterLaplacian3 = tex1Dfetch(texIn, i); 
						AverLaplacian3 = aver3D(texIn, x, y, z);
						GradLaplacian3 = grad3D(texIn, x, y, z);
					} 
					else 
					{
						CenterLaplacian3 = tex1Dfetch(texOut, i); 
						AverLaplacian3 = aver3D(texOut, x, y, z);
						GradLaplacian3 = grad3D(texOut, x, y, z);
					}
				
					if (ruidoSwitch) 
					{
						tauL = tau; 
					} 
					else 
					{
						tauL = tauC;
					}
				
					Laplaciano = AverLaplacian3 - CenterLaplacian3;
			
					Energy[i] =     	c_1 * powf(Laplaciano,2) +
										c_2 * powf(GradLaplacian3,2) + 
										c_3 * powf(CenterLaplacian3,2) + //este hay q sacar
										c_4 * powf(CenterLaplacian3,3) ;//+
										c_5 * powf(CenterLaplacian3,4);		
				}
		}
	}
}


long calc_pot2()//(float escala)
{
	int Xv, Yv, Zv;
	Xv = 0; Yv = 0; Zv = 0;

 	for (int i = 0; i < Ntotal; i++) 
	{
		
		float c_1 = hAmpl/2.0f;
		
		hPot2[i] = c_1 * (tanh(-((float)Zv - hLb) * hWd) -
					      tanh(-((float)Zv - hLa) * hWd)) + hAmpl;

		//	if(Zv <= hNz / 2)
		//{   
		//	float pot = powf(sinf((float)Zv  / escala ), 2);
		//	hPot2[i] += c_1 * ( fmin(pot, 0.4f) - 0.4f );
		//}

		/*					
		hPot2[i] = hAmpl/2.0f*(tanh(-(((float)Zv-1.0f)-hLb)*hWd) -
							tanh(-(((float)Zv-1.0f)-hLa)*hWd));
		*/									
		Xv = Xv + 1;
		if(Xv == hNx) 
		{
			Xv = 0;
			Yv = Yv+1;
			if(Yv == hNy) 
			{
				Yv = 0;
				Zv = Zv+1;
			}
		}
	
	}
	
	//save_plt(hPot2,30);
	
 return(0);
}


long calc_pot3()//(float escala)
{
	int Xv, Yv, Zv;
	Xv = 0; Yv = 0; Zv = 0;

 	for (int i = 0; i < Ntotal; i++) 
	{
		
		float c_1 = - hAmpl/2.0f;
		
		hPot3[i] = (c_1/10) * (tanh(-((float)Zv - hLb) * hWd) -
					      tanh(-((float)Zv - hLa) * hWd)) + hAmpl;

	//	if(Zv <= hNz / 2)
	//{   
	//	float pot = powf(sinf((float)Zv  / escala ), 2);
	//	hPot2[i] += c_1 * ( fmin(pot, 0.4f) - 0.4f );
	//}

		/*					
		hPot2[i] = hAmpl/2.0f*(tanh(-(((float)Zv-1.0f)-hLb)*hWd) -
							tanh(-(((float)Zv-1.0f)-hLa)*hWd));
		*/									
		Xv = Xv + 1;
		if(Xv == hNx) 
		{
			Xv = 0;
			Yv = Yv+1;
			if(Yv == hNy) 
			{
				Yv = 0;
				Zv = Zv+1;
			}
		}
	
	}
	
	//save_plt(hPot2,30);
	
 return(0);
}

//////// making geometry of system
void make_geo()
{

	// for (int i=0; i<=N-1; i++) {
		// for (int j=0; j<=N-1; j++) {
			// for (int k=0; j<=N-1; k++) {
		
		  	// X[i]=(i+1)*hL/N;
			// Y[i]=(j+1)*hL/N;
			// Z[i]=hAmp*cos(PI*nx*(i+1)/N)*cos(PI*my*(j+1)/N);
			
			// }
	
		// }
	// } 
  return;
}

//************************************************************************************************************
//									MEMORY MANAGMENT
//*************************************************************************************************************

// allocate memory for the needed arrays
long init_mem()
{
 Fini = (float *) calloc((size_t) Ntotal, sizeof(float) );
 cFinif = (float *) calloc((size_t) Ntotal, sizeof(float) );
 Energia = (float *) calloc((size_t) Ntotal, sizeof(float) );
 X = (float *) calloc((size_t) Ntotal, sizeof(float) );
 Y = (float *) calloc((size_t) Ntotal, sizeof(float) );
 Z = (float *) calloc((size_t) Ntotal, sizeof(float) );
 hPot2 = (float *) calloc((size_t) Ntotal, sizeof(float) );
 hPot3 = (float *) calloc((size_t) Ntotal, sizeof(float) );

 // ==============================================================
 if(Fini == NULL || cFinif==NULL || X==NULL || Y==NULL || Z==NULL || hPot2==NULL || Energia==NULL)
 {
	 printf("ERROR while allocating arrays \n");
	 waitKey();
	 exit(-1);
 }

 return 0;
}

//Handle the errors
static void HandleError( cudaError_t err, const char *file, int line ) 
{
    if (err != cudaSuccess) 
	{
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
		waitKey();
        exit( err );
    }
}

// clean up memory allocated on the GPU
void anim_exit( DataBlock *d ) 
{
	printf( "Ending program... \n"); 
    cudaUnbindTexture( texIn );
	cudaUnbindTexture( texPot2 );
	cudaUnbindTexture( texPot3 );
	cudaUnbindTexture( texTherF3D );
    cudaUnbindTexture( texOut );
	cudaUnbindTexture( texE );
	
    HandleErrorWrapper( cudaFree( d->dev_inSrc ) );
    HandleErrorWrapper( cudaFree( d->dev_outSrc ) );
	HandleErrorWrapper( cudaFree( d->dev_Energy ) );
	HandleErrorWrapper( cudaFree( d->Pot2 ) );
	HandleErrorWrapper( cudaFree( d->TherF3D ) );

    HandleErrorWrapper( cudaEventDestroy( d->start ) );
    HandleErrorWrapper( cudaEventDestroy( d->stop ) );
	printf( "End \n"); 
}

void run_gl( DataBlock *d) {
	dim3 blocks(hNx/Nthread, hNy/Nthread, hNz);
	dim3 threads(Nthread, Nthread);
	
	// Setting up cuRand in GPU
	if (!(etaR==0.0f)) 
	{
		printf( "Including noise..."); 
		//printf( "Setting up cuRand in GPU... "); 
		setup_kernel<<<blocks,threads>>>( d->devStates, idum ); 
		printf( "OK \n");
	}
	else
	{
		printf( "Noise is NOT included \n");
		//printf( "ERROR setting up cuRand in GPU... "); 
	}
	
	for (int k=0; k<hNF; k++) 
	{
		HandleErrorWrapper( cudaEventRecord( d->start, 0 ) );
		
		// since tex is global and bound, we have to use ha flag to
		// select which is in/out per iteration
		volatile bool dstOut = true;
		volatile bool ruidoSwitch = true;
	
		float   *out, *PhiCuadrado, *FinalPhi2, *E, *EnergiaTotal;
	
		EnergiaTotal = (float *) calloc((size_t) hNz, sizeof(float) );
		PhiCuadrado = (float *) calloc((size_t) hNz, sizeof(float) );
		FinalPhi2 = (float *) calloc((size_t) hNz, sizeof(float) );
	
		int q = d->frames+Vini+1;
	
		float FinalPhi;
		float FinalEnergy;
		
		printf( "Running... \n"); 
		for (int i=0; i<Nsteps; i++) 
		{ 
			if (dstOut) 
			{
				out = d->dev_outSrc;
			} 
			else 
			{
				out = d->dev_inSrc;
			}

			HandleErrorWrapper( cudaMemcpy( Fini, d->TherF3D, bytes, cudaMemcpyDeviceToHost ) );

			gl_kernel1<<<blocks,threads>>>(d->TherF3D, ruidoSwitch, dstOut);

			HandleErrorWrapper( cudaMemcpy( Fini, d->TherF3D, bytes, cudaMemcpyDeviceToHost ) );

			gl_kernel2<<<blocks,threads>>>(out ,d->devStates, ruidoSwitch, dstOut);

			HandleErrorWrapper( cudaMemcpy( Fini, out, bytes, cudaMemcpyDeviceToHost ) );

			dstOut = !dstOut;

			//************************************************************************************************************
			//									FOR N STEPS
			//*************************************************************************************************************
			int Nsave = fmax(1.0, Nsteps / 1000);
			if (i % Nsave == 0)
			{	
				HandleErrorWrapper( cudaMemcpy( cFinif, d->dev_inSrc, bytes, cudaMemcpyDeviceToHost ) );
				
				E = d->dev_Energy;
			
				gl_kernel4<<<blocks,threads>>>(E, ruidoSwitch, !dstOut);
			
				HandleErrorWrapper( cudaMemcpy( Energia, d->dev_Energy, bytes, cudaMemcpyDeviceToHost ) );	
					
				for(int p = 0; p < hNz; p++)
				{
					PhiCuadrado[p] = 0;
					FinalPhi2[p] = 0;
					EnergiaTotal[p] = 0;
				}
				
				float FinalPhi3  = 0.0f;
				float EnergiaTotal2 = 0.0f;
				int hNxy = hNx * hNy;
				
				for(int j = 0; j < hNz; j++)
				{		
					for(int l = 0; l < hNxy; l++)
					{				
						PhiCuadrado[j] = PhiCuadrado[j] + cFinif[ j * hNxy + l ] * cFinif[ j * hNxy + l ]; //PHI^2 BY LAYERS
						
						EnergiaTotal[j] = EnergiaTotal[j] + Energia[ j * hNxy + l ]; //saving energy by layers
					}
					
					FinalPhi2[j] = sqrtf(PhiCuadrado[j] / hNxy); //saving phi by layers
 
					FinalPhi3 = FinalPhi3 + FinalPhi2[j];//integral of phi
										
					EnergiaTotal2 = EnergiaTotal2 + (EnergiaTotal[j] / hNxy);//integral of energy
				}
			 
				FinalPhi = 0.0f;
				FinalEnergy = 0.0f;	
							
				FinalPhi = (FinalPhi3 / hNz);//Normalizing
				FinalEnergy = (EnergiaTotal2 / hNz);//Normalizing
						
				int time = 0;
				
				time = i + k * Nsteps;	
								
				if((io_error = save_dat_Energy2(FILE_PHI_TOTAL_DAT_BASENAME, FinalEnergy, FinalPhi, time, q)) != 0)
				{
					print_error("Error: saving .dat file", io_error);
				}
			}
						
			//************************************************************************************************************
			//									END BY STEPS
			//*************************************************************************************************************
		}
		
		if (etaR!=0.0f && NstepsClean>0) 
		{
			ruidoSwitch = false;
			printf( "Cleaning noise... ");
			for (int i = 0; i < NstepsClean; i++) 
			{ 
				if (dstOut) 
				{
					out = d->dev_outSrc;
				} 
				else 
				{
					out = d->dev_inSrc;
				}
				
				gl_kernel1<<<blocks,threads>>>(d->TherF3D,ruidoSwitch,dstOut);
				gl_kernel2<<<blocks,threads>>>(out,d->devStates,ruidoSwitch,dstOut);
				dstOut = !dstOut;
			}

		printf( "OK \n ");
		}
		
		//Data: device to host
		HandleErrorWrapper( cudaMemcpy( cFinif, d->dev_inSrc, bytes, cudaMemcpyDeviceToHost ) );	
		HandleErrorWrapper( cudaMemcpy( Energia, d->dev_Energy, bytes, cudaMemcpyDeviceToHost ) );				
		HandleErrorWrapper( cudaEventRecord( d->stop, 0 ) );
		HandleErrorWrapper( cudaEventSynchronize( d->stop ) );
		
		float   elapsedTime;
		HandleErrorWrapper( cudaEventElapsedTime( &elapsedTime, d->start, d->stop ) );
		d->totalTime += elapsedTime;
		++d->frames;

		printf( "File %d of %d ...",q , hNF+Vini);
		printf( "Time: %7.4f seg.\n", d->totalTime/(1000*d->frames)  );
		

		//Data saving
		printf("Saving data to filesystem. Step: %i\n", q);
		if (saveDAT==1) 
		{
			printf("Saving %s_z=%ld_%i.dat...", FILE_DAT_BASENAME, hNz, q);
			if((io_error = save_dat(FILE_DAT_BASENAME, cFinif, q)) != 0)
			{
				print_error("Error: saving .dat file", io_error);
			}
			else
			{
				printf( "OK \n");
			}
			
			if (SaveEnergy==1)
			{
				printf("Saving ENERGY %s_z=%ld_%i.dat...", FILE_E_BASENAME, hNz, q);
				if((io_error = save_dat(FILE_E_BASENAME, Energia, q)) != 0)
				{
					print_error("Error: saving .dat file", io_error);
				}
				else
				{
					printf( "OK \n");
				}
			}
		}
			
		if(k == hNF-1)
		{
			printf("Saving E TOTAL %s_z=%ld_%i.dat...", FILE_PHI_TOTAL_DAT_BASENAME, hNz, q);
			if((io_error = save_dat_Energy(FILE_PHI_TOTAL_DAT_BASENAME, FinalEnergy, FinalPhi, q)) != 0)
			{
				print_error("Error: saving .dat file", io_error);
			}
			else
			{
				printf( "OK \n");
			}
			
			printf("Saving ENERGY LAYERS %s_z=%ld_%i.dat...", FILE_PHI_DAT_BASENAME_CAPAS, hNz, q);
			if((io_error = save_dat_capas(FILE_PHI_DAT_BASENAME_CAPAS, EnergiaTotal, FinalPhi2, q)) != 0)
			{
				print_error("Error: saving .dat file", io_error);
			}
			else
			{
				printf( "OK \n");
			}
				
		}
		
		if (savePLT==1) 
		{
			printf("Saving %s_z=%ld_%i.plt...", FILE_PLT_BASENAME, hNz, q);
			if((io_error = save_plt(FILE_PLT_BASENAME, cFinif, q)) != 0)
			print_error("Error: saving .plt file", io_error);
			else
			{
				printf( "OK \n");
			}
			
			if (SaveEnergy==1)
			{
				printf("Saving ENERGY %s_z=%ld_%i.plt...", FILE_E_BASENAME, hNz, q);
				if((io_error = save_plt(FILE_E_BASENAME, Energia, q)) != 0)
				print_error("Error: saving .plt file", io_error);
				else
				{
					printf( "OK \n");
				}		
			}
		}
		
		printf("Saving data end. Step: %i\n", q);

	}

}


//************************************************************************************************************
//									FILE MANAGMENT
//*************************************************************************************************************

// GET INITIAL DATA
errno_t get_param(const char *fname)
{ 
	// read main parameters from the input file "input3D.dat"
 FILE *fp;
 char buf[100];

 errno_t error = fopen_s(&fp, fname,"r");

 if(error != 0)
	 return error;
 
  char *token;
  char  delims[] = " ,\t\n";
  char* context	 = NULL;
    
 // read all input parameters
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ",&hNx);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ",&hNy);
 
 //Read the hNZs
 //First the count
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%d ",&hNz_count);
 //If there is any hNz then read all
 if(hNz_count > 0)
 {
	hNzs = (long*) malloc(hNz_count*sizeof(long));
	for(int i = 0; i < hNz_count; i++)
	{
		fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ", &(hNzs[i]));
	}
 }
 else
 {
 return -1 ;
 }
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ",&hHold);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ",&hz_AB);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ",&z_AT);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ",&hNF);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&htau);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&htauC);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&htauP);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&ha);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&hf);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&hfP);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&hv);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&hu);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&hD);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&hDP);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&hDx);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&hDy);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&hB);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&hBP);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&etaI);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ",&etaR);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ",&idum);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ",&Nsteps);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ",&NstepsClean);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ", &hDT);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ", &hWd);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ", &hLa);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ", &hLb);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%f ", &hAmpl);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ",&hNpot);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ",&Vini);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%ld ",&Metodo);

 fgets(buf,100,fp); 
 token = strtok_s(buf, delims, &context); 
 sscanf(token,"%s", ext);

 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%d ", &saveDAT);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%d ", &savePLT);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%d ", &SaveEnergy);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%d ", &Nthread);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%s", FILE_RAND_BASENAME);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%s", FILE_PLT_BASENAME);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%s", FILE_DAT_BASENAME); 
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%s", FILE_E_BASENAME);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%s", FILE_PHI_TOTAL_DAT_BASENAME);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%s", FILE_PHI_DAT_BASENAME_CAPAS);
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%s", FILE_OPEN_BASENAME); 
 fgets(buf,100,fp);  token = strtok_s(buf, delims, &context); sscanf(token,"%s", FILE_OPEN_BASENAME_2); 
 
 fclose(fp); 


 //// print to check
 printf("Nx           = %ld, Ny=%ld NF=%ld\n", hNx, hNy, hNF);
 
 if(hNz_count > 0)
 {
	for(int i = 0; i < hNz_count; i++)
	{
 		 printf("hNz %i        = %ld \n", i, (hNzs[i]));
	}
 }
 printf("Hold            = %ld \n", hHold);
 printf("Hz_align bottom = %ld \n", hz_AB);
 printf("Hz_align top    = %ld \n", z_AT);
 printf("tau             = %f \n", htau);
 printf("tauC            = %f \n", htauC);
 printf("tauP            = %f \n", htauP);
 printf("a               = %f \n", ha);
 printf("f               = %f \n", hf);
 printf("fP              = %f \n", hfP);
 printf("v               = %f \n", hv);
 printf("u               = %f \n", hu);
 printf("D               = %f \n", hD);
 printf("DP              = %f \n", hDP);
 printf("Dx              = %f \n", hDx);
 printf("Dy              = %f \n", hDy);
 printf("B               = %f \n", hB);
 printf("BP              = %f \n", hBP);
 printf("etaI            = %f \n", etaI);
 printf("etaR            = %f \n", etaR);
 printf("idum            = %ld \n", idum);
 printf("Nsteps          = %ld \n", Nsteps);
 printf("NstepsClean     = %ld \n", NstepsClean);
 printf("DT              = %f \n", hDT);
 printf("Wd              = %f \n", hWd);
 printf("La              = %f \n", hLa);
 printf("Lb              = %f \n", hLb);
 printf("Ampl            = %f \n", hAmpl);
 printf("Npot            = %ld \n", hNpot);
 printf("Vini            = %ld \n", Vini);
 printf("Method          = %ld \n", Metodo);
 printf("Save phi.dat    = %i \n", saveDAT);
 printf("Save phi.plt    = %i \n", savePLT);
 printf("Save E.dat.plt  = %i \n", SaveEnergy); 
 printf("Extension       = %s \n", ext);
 printf("Nthread         = %d \n", Nthread);
 printf("RandName        = %s  \n", FILE_RAND_BASENAME);
 printf("PLT Name        = %s  \n", FILE_PLT_BASENAME);
 printf("DAT Name        = %s  \n", FILE_DAT_BASENAME);
 printf("Energy Name     = %s  \n", FILE_E_BASENAME);
 printf("<Phi><E>vs t    = %s  \n", FILE_PHI_TOTAL_DAT_BASENAME);
 printf("<Phi><E>vs layer= %s  \n", FILE_PHI_DAT_BASENAME_CAPAS); 
 printf("OPEN Name 1     = %s  \n", FILE_OPEN_BASENAME);
 printf("OPEN Name 2     = %s  \n", FILE_OPEN_BASENAME_2);
 fclose(fp);  
 return(0);
}

//OPEN A FILE .DAT
errno_t get_dat(char *fname)
{
	FILE *fp;

	errno_t error = fopen_s(&fp, fname,"r");
	if(error != 0)
		return error;
 
	for (int i=0; i < Ntotal; i++) 
	{
		fscanf(fp, "%f", &Fini[i]);
	}

	fclose(fp);
	return 0;
}

//OPEN A FILE .PLT
errno_t get_plt(const char *fname)
{
 FILE *fp;
 char buf[100];
 double val0, val1, val2;

 errno_t error = fopen_s(&fp, fname,"r");
 if(error != 0) 
	 return error;
 
 //Dismissing the 3 first lines
 fgets(buf,100,fp);
 fgets(buf,100,fp);
 fgets(buf,100,fp);

 	for (int i = 0; i < Ntotal; i++) {
		fscanf(fp, "%f %f %f %f", &val0, &val1, &val2, &Fini[i]); 
		//printf("Fini[%d]=%f\n",i,Fini[i]); if (i%10==0) {getchar();}
	}

 fclose(fp);

 return 0;
}

errno_t get_plt_Alineado_BOTTOM(int hz_AB, const char *fname)
{
 FILE *fp;
 char buf[100];
 double val0, val1, val2;

 int offset = hz_AB * hNx *hNy;

 errno_t error = fopen_s(&fp, fname,"r");
 if(error != 0) 
	 return error;
 
 //Dismissing the 3 first lines
 fgets(buf,100,fp);
 fgets(buf,100,fp);
 fgets(buf,100,fp);

	if (hNpot == 1) 
	{   int ipot = hLb + (hWd/2);
		for (int i = ipot ; i < (offset + ipot); i++)
		//for (int i = 0 ; i < (offset); i++) 
		{//desde comienzo del potencial
		fscanf(fp, "%f %f %f %f", &val0, &val1, &val2, &Fini[i]); 
		//printf("Fini[%d]=%f\n",i,Fini[i]); if (i%10==0) {getchar();}
		}
	}

	else
	{
		for (int i = 0 ; i < offset; i++) 
		{//desde comienzo del potencial
		fscanf(fp, "%f %f %f %f", &val0, &val1, &val2, &Fini[i]); 
		//printf("Fini[%d]=%f\n",i,Fini[i]); if (i%10==0) {getchar();}
		}
	}

 fclose(fp);

 return 0;
}

errno_t get_plt_Alineado_TOP(int hz_AB, const char *fname)
{
 FILE *fp;
 char buf[100];
 double val0, val1, val2;

 int offset = hz_AB * hNx *hNy;
 
 errno_t error = fopen_s(&fp, fname,"r");
 if(error != 0) 
	 return error;
 
 //Dismissing the 3 first lines
 fgets(buf,100,fp);
 fgets(buf,100,fp);
 fgets(buf,100,fp);

if (hNpot == 1) 
	{   int ipot = hLb + (hWd/2);
		for (int i = offset + ipot; i < Ntotal; i++) 
		//for (int i = offset; i < Ntotal; i++) 
		{
 		fscanf(fp, "%f %f %f %f", &val0, &val1, &val2, &Fini[i]); 
		}
	}

 else
	{
	 for (int i = offset; i < Ntotal; i++) 
	 {
 		fscanf(fp, "%f %f %f %f", &val0, &val1, &val2, &Fini[i]); 
 		//printf("Fini[%d]=%f\n",i,Fini[i]); if (i%10==0) {getchar();}
	}
	}

 fclose(fp);

 return 0;
}

//SAVE A .DAT
errno_t save_dat(const char *basename, float *cFini, int q)
{	FILE *fp;
	char filename[FILENAME_MAX];
	char cero[10];
	
	if (q < 10) 
		sprintf(cero, "0");
	else 
		sprintf(cero, "");
		
	sprintf_s(filename, "%s_z=%ld_%s%d.dat", basename, hNz, cero, q);
 
	errno_t error = fopen_s(&fp, filename,"w+");
	if(error != 0)
		return error;
 
	 	for (int i = 0; i < Ntotal; i++) 
		{
			fprintf(fp, "%11.6f\n", cFini[i]);
		}

	fclose(fp);

	return 0;
}


// SAVE A .PLT
errno_t save_plt(const char *basename, float *cFini, int q)
{
 FILE *fp;
 char filename[FILENAME_MAX];
 char cero[10];
	
	if (q < 10) 
		sprintf(cero, "0");
	else 
		sprintf(cero, "");
		
 sprintf_s(filename, "%s_z=%ld_%s%d.plt", basename, hNz, cero, q);

 errno_t error = fopen_s(&fp, filename, "w+");
 if(error != 0)
	 return error;

 char header1[100];
 char header2[100];
 char header3[100];

 sprintf_s(header1, "TITLE =\"Pattern 2D\"");
 sprintf_s(header2, "VARIABLES = \"X\",\"Y\",\"Z\",\"Psi\"");
 sprintf_s(header3, "ZONE I=%d, J=%d, K=%d, F=Point", hNx, hNy, hNz);
 
 fprintf(fp, "%s\n", header1);
 fprintf(fp, "%s\n", header2);
 fprintf(fp, "%s\n", header3);
 
 int Nxy = hNx * hNy;
 
 for(int k = 0; k < hNz; k++){
	for(int i = 0; i < hNy; i++){
		for(int j = 0; j < hNx; j++){

			fprintf(fp, "%3d.00 %3d.00 %3d.000000 %1.6f\n", i, j, k, cFini[k * Nxy + j * hNy + i]);
		}
	}
}
	
 fclose(fp);

 return 0;
}

// SAVE A ENERGY VS ESPESOR
errno_t save_dat_Energy2(const char *basename, float EnergiaTotal, float FinalPhi, int k, int q)
{	
	FILE *fp;
	char filename[FILENAME_MAX];
	char cero[10];
	
	if (q < 10) 
		sprintf(cero, "0");
	else 
		sprintf(cero, "");
		
	sprintf_s(filename, "%s_z=%ld_all.dat", basename, hNz);
	
	errno_t error = fopen_s(&fp, filename,"a");
	
	
	if(error != 0)
		return error;
		
			fprintf(fp, "%ld %ld %11.6f %11.6f\n", hNz, k, EnergiaTotal, FinalPhi);

	fclose(fp);

	return 0;
}

// SAVE A ENERGY VS ESPESOR
errno_t save_dat_Energy(const char *basename, float EnergiaTotal, float FinalPhi, int q)
{	
	FILE *fp;
	char filename[FILENAME_MAX];
	char cero[10];
	
	if (q < 10) 
		sprintf(cero, "0");
	else 
		sprintf(cero, "");
		
	sprintf_s(filename, "%s.dat", basename);
 
		
	errno_t error = fopen_s(&fp, filename,"a");
	
	
	if(error != 0)
		return error;
		
			fprintf(fp, "%ld %11.6f %11.6f\n", hNz, EnergiaTotal, FinalPhi);

	fclose(fp);

	return 0;
}

//SAVE A .DAT ENERGY POR CAPAS
errno_t save_dat_capas(const char *basename, float *energy, float *phi, int q)
{	
	FILE *fp;
	char filename[FILENAME_MAX];
	char cero[10];
	
	if (q < 10) 
		sprintf(cero, "0");
	else 
		sprintf(cero, "");
		
	sprintf_s(filename, "%s_z=%ld_%s%d.dat", basename, hNz, cero, q);
 
	errno_t error = fopen_s(&fp, filename,"w+");
	
	if(error != 0)
		return error;
 
	float Nxy = hNx*hNy;
	for(int j = 0; j < hNz; j++)
	{
		fprintf(fp, "%d %11.6f %11.6f\n", j, energy[j]/(Nxy), phi[j]);
	}
		

	fclose(fp);

	return 0;
}


//GENERATE RANDOM MATRIX FINI 
errno_t gen_rand(const char* fname, float *Fini){
 FILE *fp;
 
	char filename[FILENAME_MAX];
	sprintf_s(filename, "%s.dat", fname);

	errno_t error = fopen_s(&fp, filename, "w+");

	if(error != 0)
		return error;
 
	for (int i=0; i < Ntotal; i++) 
	{
		Fini[i] = etaI*(1.0f - 2.0f * ran1(&idum));
		fprintf(fp, "%15.6E\n", Fini[i]);
	}
	
 fclose(fp);
 return 0;
}

errno_t gen_rand_Alineado_BOTTOM(int hz_AB, const char* fname, float *Fini){
 FILE *fp;
 
	char filename[FILENAME_MAX];
	sprintf_s(filename, "%s_alineado.dat", fname);

	errno_t error = fopen_s(&fp, filename, "w+");

	if(error != 0)
		return error; 
 
	int offset = hz_AB * hNx * hNy;
 
 	for (int i = offset; i < Ntotal; i++) 
	{
		Fini[i] = etaI*(1.0f - 2.0f * ran1(&idum));
		fprintf(fp, "%15.6E\n", Fini[i]);
	}
	
 fclose(fp);
 return 0;
}

//PRINT ERROR
void print_error(const char *msg, int err) 
{
	printf("%s .Errno %i\n", msg, err);
	printf("Press any key to exit\n");
	char buf[10];
	scanf_s(buf);
	exit(err);
}

//WAIT BEFORE EXIT
void waitKey()
{
	printf("Press any key to continue...\n");
	char buf[10];
	scanf_s(buf);
}