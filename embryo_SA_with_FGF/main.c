#include <stdio.h>
#include <assert.h>
#include <sys/sysctl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <math.h>

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

/* Parameter setup for main variants */
// Selective adhesion (SA) mechanism setup (should be a string)
#define SA_type "0" /* "0":Epha4/Efnb2 SA; "1":no SA; "2": symmetric SA;
                    "3": None-biased  asymmetric SA
                    "4": EPi-biased asymmetric SA
                    "5": PE-biased asymmetric SA*/
// Selective adhesion (SA) on time
#define SA_on_time 95000 /* start from 32C: 75000
                            start from 64C: 95000
                            start from 128C: 120000 */
// FGF on time
#define FGF_on_time 0 /* on from beginning: 0
                        on from 16C: 39000
                        on from 32C: 75000
                        on from 64C: 95000
                        on from 128C: 120000 */
// FGF off time
#define FGF_off_time 160000 /* off at 32C: 75000
                                off at 64C: 95000
                                off at late 128C: 160000
                                never off: 260000 */
// output file name (space file end with "XXX_space.timestamp.PQR" in "PQR" folder;
//                   Gene file end with "XXX_gene.timestamp.txt" in "signal" folder)
#define File_name "example_file"
/* End of parameter setup*/

#define PI 3.1415926535
#define BC 0
#define RATIO 1
#define NX 50
#define NY NX
#define NZ NX/RATIO
#define NEQ  NX*NY*NZ          /* problem dimension */

#define NO 0
#define YES 1
#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

#define LX 60.0
#define LY LX
#define LZ LX/RATIO
#define DX LX/NX
#define DY DX
#define DZ DX

#define GROWTH_FR 1000
#define DIVISION_FR 10000
#define OUTPUT_FR 1000

#define MAX_ITERATIONS 260001
#define MAX_CHEM_ITER  1000
#define MAX_CELL 1280
#define MAX_ELE 20

static const int nx = NX;
static const int ny = NY;
static const int nz = NZ;
float lx = 13.00f;
static const float ly = LY;
static const float lz = LZ;

extern int errno;

/* Private functions */
float ***chem1, ***chem2, ***chem_diff;
float **signal;
int *id_is_out;
int *id_is_f;

int ***grid, ***grid1, ***grid2, ***grid3;
void  initial_chem_paras(float * x, float * y, float * z, int * id, int * cell_type, int * ele_type, int cell_no, int * ele_per_cell);
void  cell_in_block(float * x, float * y, float * z, int * id, int * cell_type, int * ele_type, int cell_no, int * ele_per_cell, int * judgeout, int * judgein);

void  signal_init();
void  signal_update(int, float, int, int*, int*, int, int*, int*, int*, int*, int*);


float min(float a, float b)
{
    return (a > b)? b : a;
}
float max(float a, float b)
{
    return (a > b)? a : b;
}


float RNG()
{
    float random;
    random=(1.0*(rand()%RAND_MAX))/(1.0*RAND_MAX);
    return (random);
}

float RNGn(float mu, float sigma)
{
  float U1, U2, W, mult;
  static float X1, X2;
  static int call = 0;

  if (call == 1)
  {
      call = !call;
      return (mu + sigma * (float) X2);
  }

  do
  {
      U1 = -1 + ((float) rand () / RAND_MAX) * 2;
      U2 = -1 + ((float) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
  }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (float) X1);

}

float distance3(float x1, float y1, float z1, float x2, float y2, float z2)
{
    float dx = fabs(x1 - x2);
    float dy = fabs(y1 - y2);
    float dz = fabs(z1 - z2);
    
    return sqrt(dx*dx + dy*dy + dz*dz) + 0.0001;
}

void InitialReader(const char InName[], float * x, float * y, float * z, int * id, int * cell_type, int * ele_type, int * cell_no, int * ele_no, int * ele_per_cell)
{
    FILE *parIn;
    int i, j;
    float px, py, pz;
    int pid, ptype, tmp;
    int flag_return;

    parIn = fopen(InName , "r");

    if (parIn == NULL)
    {
       printf("Parameter File %s was not open\n", InName);
       exit(1);
    }
    printf("\nEnter read_initial_conditions():\n");
    printf("filename=%s\n", InName);

    flag_return = fscanf(parIn, "%d", cell_no);
    printf("cell_no = %d\n", *cell_no);
    flag_return = fscanf(parIn, "%d", ele_no);
    printf("ele_no = %d\n", *ele_no);
    tmp = (*ele_no)/(*cell_no);
    (*ele_no) = 0;
    for(i = 0; i < *cell_no; i++)
    {
    for(j = i*MAX_ELE+0; j < i*MAX_ELE+tmp; j++)
    {
        flag_return = fscanf(parIn, "%d %f %f %f %d", &pid, &px, &py, &pz, &ptype);
        printf("%d %f %f %f %d\n", pid, px, py, pz, ptype);
        x[j] = px;
        y[j] = py;
        z[j] = pz;
        id[j] = pid;
        ele_per_cell[pid-1]++;
        ele_type[j] = ptype;
        (*ele_no)++;
    }
    cell_type[i] = 1;
    }

    fclose(parIn);
}

void output(float * x, float * y, float * z, int * id, int * cell_type, int * ele_type, int cell_no, int ele_no, int * ele_per_cell, int step, float * chem1_gpu, float * chem2_gpu, float * xc, float * yc, float *zc, int * flag_outer, int * id_is_out, float ** signal, int * colevel, int * id_is_f)
{
    int i, k, check;
    char filename[50], outdir[50];
    int cellt;
    FILE *fp = NULL;

    struct stat st;

    /////////////////////////////////////////////////////////
    // PQR
    /////////////////////////////////////////////////////////
    sprintf(outdir, "%s", "PQR");

    if(stat(outdir, &st) == -1)
    {
        if(errno == ENOENT)
        {
        check = mkdir(outdir, S_IRWXU);
        if(check != 0)
        {
        (void) printf("WARNING in output(), directory "
                     "%s doesn't exist and can't be created\n",outdir);
        }
        else
        {
            printf("created the dirctory %s\n", outdir);
        }
        }
    }
    sprintf(filename, "%s/" File_name "_space.%d.pqr", outdir, step);
    fp = fopen(filename, "w");
    if(fp == NULL)
    {
        printf("Failed to open file for writing\n");
        exit(-1);
    }
    for(i=0;i<ele_no;i++)
    {
        k = i;
        cellt=id_is_out[id[k]-1]+1;
        fprintf(fp, "ATOM %6d C  THR   %5d      % .3f % .3f % .3f % .4f % .4f\n",
                k+1, id[k]+cellt*1000+id_is_f[id[k]-1]*200, x[k], y[k], z[k], 0.1*id_is_out[id[k]-1], 0.95);
    }
    fclose(fp);
    
    sprintf(outdir, "%s", "CO");
    
    if(stat(outdir, &st) == -1)
    {
        if(errno == ENOENT)
        {
            check = mkdir(outdir, S_IRWXU);
            if(check != 0)
            {
                (void) printf("WARNING in output(), directory "
                              "%s doesn't exist and can't be created\n",outdir);
            }
            else
            {
                printf("created the dirctory %s\n", outdir);
            }
        }
    }
    
    
    sprintf(filename, "%s/" File_name "_gene.%d.pqr", outdir, step);
    fp = fopen(filename, "w");
    if(fp == NULL)
    {
        printf("Failed to open file for writing\n");
        exit(-1);
    }
    for(i=0;i<ele_no;i++)
    {
        k = i;
        cellt=colevel[id[k]-1]+1;
        fprintf(fp, "ATOM %6d C  THR   %5d      % .3f % .3f % .3f % .4f % .4f\n",
                k+1, id[k]+cellt*1000, x[k], y[k], z[k], 0.1*id_is_out[id[k]-1], 0.95);
    }
    fclose(fp);
    
    sprintf(outdir, "%s", "signal");
    
    if(stat(outdir, &st) == -1)
    {
        if(errno == ENOENT)
        {
            check = mkdir(outdir, S_IRWXU);
            if(check != 0)
            {
                (void) printf("WARNING in output(), directory "
                              "%s doesn't exist and can't be created\n",outdir);
            }
            else
            {
                printf("created the dirctory %s\n", outdir);
            }
        }
    }
    
    sprintf(filename, "%s/" File_name "_gene.%d.txt", outdir, step);
    fp = fopen(filename, "w");
    if(fp == NULL)
    {
        printf("Failed to open file for writing\n");
        exit(-1);
    }
    for(i=0;i<cell_no;i++)
    {
        k = i;
        fprintf(fp, "%6d %6d % .3f % .3f % .3f % .3f % .3f % .3f % .3f % .3f\n",
                k+1, id_is_out[i], signal[i][0], signal[i][1], signal[i][2], signal[i][3], signal[i][4], signal[i][5], signal[i][6], signal[i][7]);
    }
    fclose(fp);


    if(step == 260000)
    {
    sprintf(outdir, "%s", "signal");
    
    if(stat(outdir, &st) == -1)
    {
        if(errno == ENOENT)
        {
            check = mkdir(outdir, S_IRWXU);
            if(check != 0)
            {
                (void) printf("WARNING in output(), directory "
                              "%s doesn't exist and can't be created\n",outdir);
            }
            else
            {
                printf("created the dirctory %s\n", outdir);
            }
        }
    }
    
    sprintf(filename, "%s/initialsignal.txt", outdir);
    fp = fopen(filename, "w");
    if(fp == NULL)
    {
        printf("Failed to open file for writing\n");
        exit(-1);
    }
    for(i=0;i<cell_no;i++)
    {
        k = i;
        fprintf(fp, "Cell %6d with chemical  % .3f % .3f % .3f % .3f\n",
                k+1, signal[i][0], signal[i][1], signal[i][2], signal[i][3]);
    }
    fclose(fp);
    }

}

#pragma mark -
#pragma mark Utilities
char * load_program_source(const char *filename)
{
	
	struct stat statbuf;
	FILE *fh;
	char *source;
	int flag_return;

	fh = fopen(filename, "r");
	if (fh == 0)
		return 0; 
	
	stat(filename, &statbuf);
	source = (char *) malloc(statbuf.st_size + 1);
	flag_return = fread(source, statbuf.st_size, 1, fh);
	source[statbuf.st_size] = '\0'; 
	
	return source; 
} 

#pragma mark -
#pragma mark Main OpenCL Routine
int runCL(float * x, float * y, float * z, int * id, int * cell_type, int * ele_type, float dt, 
          int max_ele_no, 
          int ele_no, int cell_no, int * ele_per_cell)
{
    int flag_print = 0;

	cl_program program[10];
	cl_kernel kernel[10];
	
	cl_command_queue cmd_queue;
	cl_context   context;
	
	cl_int err = 0;
	size_t returned_size = 0;
	size_t buffer_size;
	
	// Allocate some memory and a place for the results
	float * xF = (float *)malloc(max_ele_no*sizeof(float));
	float * yF = (float *)malloc(max_ele_no*sizeof(float));
	float * zF = (float *)malloc(max_ele_no*sizeof(float));
	float * xc = (float *)malloc(MAX_CELL*sizeof(float));
	float * yc = (float *)malloc(MAX_CELL*sizeof(float));
	float * zc = (float *)malloc(MAX_CELL*sizeof(float));
	int * flag = (int *)malloc(MAX_CELL*sizeof(int));
    float * chem1_gpu = (float *)malloc(NX*NY*NZ*sizeof(float));
    float * chem2_gpu = (float *)malloc(NX*NY*NZ*sizeof(float));
    float * chem11_gpu = (float *)malloc(NX*NY*NZ*sizeof(float));
    float * chem22_gpu = (float *)malloc(NX*NY*NZ*sizeof(float));
    int * grid_gpu = (int *)malloc(NX*NY*NZ*sizeof(int));
    int * grid1_gpu = (int *)malloc(NX*NY*NZ*sizeof(int));
    int * grid2_gpu = (int *)malloc(NX*NY*NZ*sizeof(int));
    int * grid3_gpu = (int *)malloc(NX*NY*NZ*sizeof(int));
    int * flag_outer = (int *)malloc(MAX_CELL*sizeof(int));
    int * flag_contact = (int *)malloc(MAX_CELL*sizeof(int));
    float * cycle = (float *)malloc(MAX_CELL*sizeof(float));
    float * randx = (float *)malloc(max_ele_no*sizeof(float));
    float * randy = (float *)malloc(max_ele_no*sizeof(float));
    float * randz = (float *)malloc(max_ele_no*sizeof(float));
    int * change_type = (int *)malloc(MAX_CELL*sizeof(int));
    float * change_time = (float *)malloc(MAX_CELL*sizeof(float));
    int * judgeout = (int *)malloc(MAX_CELL*sizeof(int));
    int * judgein = (int *)malloc(MAX_CELL*sizeof(int));
    int * colevel = (int *)malloc(MAX_CELL*sizeof(int));
    int * mother_id = (int *)malloc(MAX_CELL*sizeof(int));


    srand(time(NULL));

    for(int i = 0; i < max_ele_no; i++)
    {
        xF[i] = yF[i] = zF[i] = 0.0f;
        randx[i] = randy[i] = randz[i] =0.0f;
    }
    for(int i = 0; i < MAX_CELL; i++)
    {
        xc[i] = yc[i] = zc[i] = 0.0;
        flag[i] = flag_outer[i] = 0;
        cycle[i] = 0;
        change_type[i]=0;
        change_time[i]=0;
        judgeout[i]=50;
        judgein[i]=0;
        colevel[i]=0;
        mother_id[i]=2;
    }

    unsigned int * rng_x = (unsigned int *)malloc(max_ele_no*sizeof(unsigned int));
    unsigned int * rng_c = (unsigned int *)malloc(max_ele_no*sizeof(unsigned int));
    for(int i = 0; i < max_ele_no; i++)
    {
        rng_x[i] = rand();
        rng_c[i] = rand();
    }

	cl_mem x_mem, y_mem, z_mem;
	cl_mem xc_mem, yc_mem, zc_mem;
    cl_mem id_mem, cell_type_mem, ele_type_mem;
	cl_mem xF_mem, yF_mem, zF_mem;
    cl_mem randx_mem, randy_mem, randz_mem;
    cl_mem ele_per_cell_mem;
    cl_mem flag_mem, flag_outer_mem, flag_contact_mem, id_is_out_mem, colevel_mem;
    cl_mem chem1_gpu_mem, chem2_gpu_mem;
    cl_mem chem11_gpu_mem, chem22_gpu_mem;
    cl_mem grid_gpu_mem, grid1_gpu_mem, grid2_gpu_mem, grid3_gpu_mem;
    cl_mem rng_x_mem, rng_c_mem;

    // Get platform and device information
    cl_platform_id *platform_id = NULL;
    cl_device_id cpu = NULL, device = NULL;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;

    // get all platforms
    platform_id = (cl_platform_id*) malloc(sizeof(cl_platform_id) * 2);
    clGetPlatformIDs(2, platform_id, NULL);
    err = clGetDeviceIDs( platform_id[0], CL_DEVICE_TYPE_CPU, 1,
            &device, &ret_num_devices);
    printf("err = %d\n", err);
    assert(err == CL_SUCCESS);


#pragma mark Device Information
	{
		cl_char vendor_name[1024] = {0};
		cl_char device_name[1024] = {0};
		cl_char buf[1024] = {0};
		err = clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(vendor_name), 
							  vendor_name, &returned_size);
		err |= clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_name), 
							  device_name, &returned_size);
        err |= clGetDeviceInfo(device, CL_DEVICE_VERSION, sizeof(buf), 
							  buf, &returned_size);
		assert(err == CL_SUCCESS);
		printf("Connecting to %s %s supporting ", vendor_name, device_name);
        printf("%s...\n", buf);
	}
	
#pragma mark Context and Command Queue
	{
		context = clCreateContext(0, 1, &device, NULL, NULL, &err);
		assert(err == CL_SUCCESS);
		
		cmd_queue = clCreateCommandQueue(context, device, 0, NULL);
	}
	
#pragma mark Program and Kernel Creation
	{
		const char * filename = "movement_" SA_type ".cl";
		char *program_source = load_program_source(filename);

		program[0] = clCreateProgramWithSource(context, 1, (const char**)&program_source,
											   NULL, &err);
		assert(err == CL_SUCCESS);

		err = clBuildProgram(program[0], 0, NULL, NULL, NULL, NULL);
		assert(err == CL_SUCCESS);

        char build[2048];
        clGetProgramBuildInfo(program[0], device, CL_PROGRAM_BUILD_LOG, 2048, build, NULL);
        printf("build log=%s\n", build);	

		kernel[0] = clCreateKernel(program[0], "movement", &err);

        size_t workgroup_size;
        err = clGetKernelWorkGroupInfo(kernel[0], device, CL_KERNEL_WORK_GROUP_SIZE,
                                sizeof(size_t), &workgroup_size, NULL);
        printf("workgroup_size=%zd\n",workgroup_size);
	}
    
    {
        const char * filename = "cellcontact.cl";
        char *program_source = load_program_source(filename);
        
        program[2] = clCreateProgramWithSource(context, 1, (const char**)&program_source,
                                               NULL, &err);
        assert(err == CL_SUCCESS);
        
        err = clBuildProgram(program[2], 0, NULL, NULL, NULL, NULL);
        assert(err == CL_SUCCESS);
        
        char build[2048];
        clGetProgramBuildInfo(program[2], device, CL_PROGRAM_BUILD_LOG, 2048, build, NULL);
        printf("build log=%s\n", build);
        
        kernel[2] = clCreateKernel(program[2], "cellcontact", &err);
        
        size_t workgroup_size;
        err = clGetKernelWorkGroupInfo(kernel[2], device, CL_KERNEL_WORK_GROUP_SIZE,
                                       sizeof(size_t), &workgroup_size, NULL);
        printf("workgroup_size=%zd\n",workgroup_size);
    }

	
    {
        const char * filename = "cellcenter.cl";
        char *program_source = load_program_source(filename);
        // printf("program_source=\n%s\n", program_source);

        program[4] = clCreateProgramWithSource(context, 1, (const char**)&program_source,
                                               NULL, &err);
        assert(err == CL_SUCCESS);
        // printf("program=\n%s\n",program[0]);

        err = clBuildProgram(program[4], 0, NULL, NULL, NULL, NULL);
        //assert(err == CL_SUCCESS);

        char build[2048];
        clGetProgramBuildInfo(program[4], device, CL_PROGRAM_BUILD_LOG, 2048, build, NULL);
        printf("build log=%s\n", build);

        // Now create the kernel "objects" that we want to use in the example file 
        kernel[4] = clCreateKernel(program[4], "cellcenter", &err);

        size_t workgroup_size;
        err = clGetKernelWorkGroupInfo(kernel[4], device, CL_KERNEL_WORK_GROUP_SIZE,
                                sizeof(size_t), &workgroup_size, NULL);
        printf("workgroup_size=%zd\n",workgroup_size);
    }

#pragma mark Memory Allocation
	{
		// Allocate memory on the device to hold our data and store the results into
		buffer_size = sizeof(float) * max_ele_no;
		
		x_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err = clEnqueueWriteBuffer(cmd_queue, x_mem, CL_TRUE, 0, buffer_size,
								   (void*)x, 0, NULL, NULL);
		y_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err |= clEnqueueWriteBuffer(cmd_queue, y_mem, CL_TRUE, 0, buffer_size,
									(void*)y, 0, NULL, NULL);
		z_mem	= clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err |= clEnqueueWriteBuffer(cmd_queue, z_mem, CL_TRUE, 0, buffer_size,
									(void*)z, 0, NULL, NULL);
		assert(err == CL_SUCCESS);

		xF_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err = clEnqueueWriteBuffer(cmd_queue, xF_mem, CL_TRUE, 0, buffer_size,
								   (void*)xF, 0, NULL, NULL);
		yF_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err |= clEnqueueWriteBuffer(cmd_queue, yF_mem, CL_TRUE, 0, buffer_size,
									(void*)yF, 0, NULL, NULL);
		zF_mem	= clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err |= clEnqueueWriteBuffer(cmd_queue, zF_mem, CL_TRUE, 0, buffer_size,
									(void*)zF, 0, NULL, NULL);
        randx_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err = clEnqueueWriteBuffer(cmd_queue, randx_mem, CL_TRUE, 0, buffer_size,
                                   (void*)xF, 0, NULL, NULL);
        randy_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, randy_mem, CL_TRUE, 0, buffer_size,
                                    (void*)yF, 0, NULL, NULL);
        randz_mem	= clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, randz_mem, CL_TRUE, 0, buffer_size,
                                    (void*)zF, 0, NULL, NULL);
		assert(err == CL_SUCCESS);

		buffer_size = sizeof(int) * max_ele_no;
		id_mem	= clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err = clEnqueueWriteBuffer(cmd_queue, id_mem, CL_TRUE, 0, buffer_size,
									(void*)id, 0, NULL, NULL);
		ele_type_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err |= clEnqueueWriteBuffer(cmd_queue, ele_type_mem, CL_TRUE, 0, buffer_size,
									(void*)ele_type, 0, NULL, NULL);
		assert(err == CL_SUCCESS);
		
		buffer_size = sizeof(int) * MAX_CELL;
        cell_type_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err = clEnqueueWriteBuffer(cmd_queue, cell_type_mem, CL_TRUE, 0, buffer_size,
                                    (void*)cell_type, 0, NULL, NULL);
		ele_per_cell_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err |= clEnqueueWriteBuffer(cmd_queue, ele_per_cell_mem, CL_TRUE, 0, buffer_size,
									(void*)ele_per_cell, 0, NULL, NULL);
		flag_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err |= clEnqueueWriteBuffer(cmd_queue, flag_mem, CL_TRUE, 0, buffer_size,
									(void*)flag, 0, NULL, NULL);
        flag_outer_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, flag_outer_mem, CL_TRUE, 0, buffer_size,
                                    (void*)flag_outer, 0, NULL, NULL);
        flag_contact_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, flag_outer_mem, CL_TRUE, 0, buffer_size,
                                    (void*)flag_outer, 0, NULL, NULL);
        id_is_out_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, id_is_out_mem, CL_TRUE, 0, buffer_size,
                                    (void*)id_is_out, 0, NULL, NULL);
        colevel_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, colevel_mem, CL_TRUE, 0, buffer_size,
                                    (void*)colevel, 0, NULL, NULL);
		assert(err == CL_SUCCESS);

		buffer_size = sizeof(float) * MAX_CELL;
		xc_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err = clEnqueueWriteBuffer(cmd_queue, xc_mem, CL_TRUE, 0, buffer_size,
								   (void*)xc, 0, NULL, NULL);
		yc_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err |= clEnqueueWriteBuffer(cmd_queue, yc_mem, CL_TRUE, 0, buffer_size,
								   (void*)yc, 0, NULL, NULL);
		zc_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
		err |= clEnqueueWriteBuffer(cmd_queue, zc_mem, CL_TRUE, 0, buffer_size,
								   (void*)zc, 0, NULL, NULL);
		assert(err == CL_SUCCESS);

        buffer_size = sizeof(float) * NX*NY*NZ;
        chem1_gpu_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err = clEnqueueWriteBuffer(cmd_queue, chem1_gpu_mem, CL_TRUE, 0, buffer_size,
                                   (void*)chem1_gpu, 0, NULL, NULL);
        chem2_gpu_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, chem2_gpu_mem, CL_TRUE, 0, buffer_size,
                                   (void*)chem2_gpu, 0, NULL, NULL);
        chem11_gpu_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err = clEnqueueWriteBuffer(cmd_queue, chem11_gpu_mem, CL_TRUE, 0, buffer_size,
                                   (void*)chem11_gpu, 0, NULL, NULL);
        chem22_gpu_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, chem22_gpu_mem, CL_TRUE, 0, buffer_size,
                                   (void*)chem22_gpu, 0, NULL, NULL);
        buffer_size = sizeof(int) * NX*NY*NZ;
        grid_gpu_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, grid_gpu_mem, CL_TRUE, 0, buffer_size,
                                   (void*)grid_gpu, 0, NULL, NULL);
        grid1_gpu_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, grid1_gpu_mem, CL_TRUE, 0, buffer_size,
                                   (void*)grid1_gpu, 0, NULL, NULL);
        grid2_gpu_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, grid2_gpu_mem, CL_TRUE, 0, buffer_size,
                                   (void*)grid2_gpu, 0, NULL, NULL);
        grid3_gpu_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, grid3_gpu_mem, CL_TRUE, 0, buffer_size,
                                   (void*)grid3_gpu, 0, NULL, NULL);
        assert(err == CL_SUCCESS);

        buffer_size = max_ele_no*sizeof(unsigned int);
        rng_x_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err = clEnqueueWriteBuffer(cmd_queue, rng_x_mem, CL_TRUE, 0, buffer_size,
                                   (void*)rng_x, 0, NULL, NULL);
        rng_c_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, rng_c_mem, CL_TRUE, 0, buffer_size,
                                   (void*)rng_c, 0, NULL, NULL);
        assert(err == CL_SUCCESS);

		// Get all of the stuff written and allocated 
		clFinish(cmd_queue);
	}
	
#pragma mark Kernel Arguments
	{
        int max_ele = MAX_ELE;
        int i_cycle = 0;
        int flag_division = 0, id1 = 0, id2 = 0;
        float xcc = 0.0f, ycc = 0.0f, zcc = 1.0f;
        int sa_on_time = SA_on_time;
		// Now setup the arguments to kernel "movement"
		err  = clSetKernelArg(kernel[0],  0, sizeof(cl_mem), &id_mem);
		err |= clSetKernelArg(kernel[0],  1, sizeof(cl_mem), &x_mem);
		err |= clSetKernelArg(kernel[0],  2, sizeof(cl_mem), &y_mem);
		err |= clSetKernelArg(kernel[0],  3, sizeof(cl_mem), &z_mem);
		err |= clSetKernelArg(kernel[0],  4, sizeof(cl_mem), &ele_type_mem);
		err |= clSetKernelArg(kernel[0],  5, sizeof(int), &cell_no);
		err |= clSetKernelArg(kernel[0],  6, sizeof(int), &ele_no);
		err |= clSetKernelArg(kernel[0],  7, sizeof(cl_mem), &ele_per_cell_mem);
		err |= clSetKernelArg(kernel[0],  8, sizeof(cl_mem), &xc_mem);
		err |= clSetKernelArg(kernel[0],  9, sizeof(cl_mem), &yc_mem);
		err |= clSetKernelArg(kernel[0], 10, sizeof(cl_mem), &zc_mem);
		err |= clSetKernelArg(kernel[0], 11, sizeof(float), &dt);
		err |= clSetKernelArg(kernel[0], 12, sizeof(cl_mem), &xF_mem);
		err |= clSetKernelArg(kernel[0], 13, sizeof(cl_mem), &yF_mem);
		err |= clSetKernelArg(kernel[0], 14, sizeof(cl_mem), &zF_mem);
		err |= clSetKernelArg(kernel[0], 15, sizeof(float), &lx);
		err |= clSetKernelArg(kernel[0], 16, sizeof(float), &ly);
		err |= clSetKernelArg(kernel[0], 17, sizeof(float), &lz);
		err |= clSetKernelArg(kernel[0], 18, sizeof(int), &max_ele);
		err |= clSetKernelArg(kernel[0], 19, sizeof(cl_mem), &cell_type_mem);
		err |= clSetKernelArg(kernel[0], 20, sizeof(int), &i_cycle);
        err |= clSetKernelArg(kernel[0], 21, sizeof(cl_mem), &id_is_out_mem);
        err |= clSetKernelArg(kernel[0], 22, sizeof(int), &i_cycle);
        err |= clSetKernelArg(kernel[0], 23, sizeof(int), &flag_division);
        err |= clSetKernelArg(kernel[0], 24, sizeof(int), &id1);
        err |= clSetKernelArg(kernel[0], 25, sizeof(int), &id2);
        err |= clSetKernelArg(kernel[0], 26, sizeof(cl_mem), &randx_mem);
        err |= clSetKernelArg(kernel[0], 27, sizeof(cl_mem), &randy_mem);
        err |= clSetKernelArg(kernel[0], 28, sizeof(cl_mem), &randz_mem);
        err |= clSetKernelArg(kernel[0], 29, sizeof(cl_mem), &colevel_mem);
        err |= clSetKernelArg(kernel[0], 30, sizeof(float), &xcc);
        err |= clSetKernelArg(kernel[0], 31, sizeof(float), &ycc);
        err |= clSetKernelArg(kernel[0], 32, sizeof(float), &zcc);
        err |= clSetKernelArg(kernel[0], 33, sizeof(int), &sa_on_time);
		assert(err == CL_SUCCESS);
	}
    
    {
        err  = clSetKernelArg(kernel[2],  0, sizeof(cl_mem), &xc_mem);
        err |= clSetKernelArg(kernel[2],  1, sizeof(cl_mem), &yc_mem);
        err |= clSetKernelArg(kernel[2],  2, sizeof(cl_mem), &zc_mem);
        err |= clSetKernelArg(kernel[2],  3, sizeof(float), &lx);
        err |= clSetKernelArg(kernel[2],  4, sizeof(float), &ly);
        err |= clSetKernelArg(kernel[2],  5, sizeof(float), &lz);
        err |= clSetKernelArg(kernel[2],  6, sizeof(cl_mem), &flag_contact_mem);
        err |= clSetKernelArg(kernel[2],  7, sizeof(int), &cell_no);
        err |= clSetKernelArg(kernel[2],  8, sizeof(cl_mem), &id_is_out_mem);
        assert(err == CL_SUCCESS);
    }
	
    {
        int max_ele = MAX_ELE;
        int tmp_i = 0;
        // Now setup the arguments to kernel "cellcenter"
        err  = clSetKernelArg(kernel[4],  0, sizeof(cl_mem), &id_mem);
        err |= clSetKernelArg(kernel[4],  1, sizeof(cl_mem), &x_mem);
        err |= clSetKernelArg(kernel[4],  2, sizeof(cl_mem), &y_mem);
        err |= clSetKernelArg(kernel[4],  3, sizeof(cl_mem), &z_mem);
        err |= clSetKernelArg(kernel[4],  4, sizeof(cl_mem), &ele_per_cell_mem);
        err |= clSetKernelArg(kernel[4],  5, sizeof(cl_mem), &xc_mem);
        err |= clSetKernelArg(kernel[4],  6, sizeof(cl_mem), &yc_mem);
        err |= clSetKernelArg(kernel[4],  7, sizeof(cl_mem), &zc_mem);
		err |= clSetKernelArg(kernel[4],  8, sizeof(int), &max_ele);
		err |= clSetKernelArg(kernel[4],  9, sizeof(float), &lx);
		err |= clSetKernelArg(kernel[4], 10, sizeof(float), &ly);
		err |= clSetKernelArg(kernel[4], 11, sizeof(float), &lz);
		err |= clSetKernelArg(kernel[4], 12, sizeof(int), &ele_no);
        err |= clSetKernelArg(kernel[4], 13, sizeof(cl_mem), &flag_outer_mem);
        err |= clSetKernelArg(kernel[4], 14, sizeof(int), &tmp_i);
        assert(err == CL_SUCCESS);
    }

#pragma mark Execution and Read
    for(int i = 0; i < NX; i++)    
        for(int j = 0; j < NY; j++)    
            for(int k = 0; k < NZ; k++) 
    {
        int indx = i + j*NX + k*NX*NY;
        chem1_gpu[indx] = chem1[i][j][k];
        chem2_gpu[indx] = chem2[i][j][k];
        grid_gpu[indx] = grid[i][j][k];
        grid1_gpu[indx] = grid1[i][j][k];
        grid2_gpu[indx] = grid2[i][j][k];
        grid3_gpu[indx] = grid3[i][j][k];
    }   
    err = clEnqueueWriteBuffer(cmd_queue, chem1_gpu_mem, CL_TRUE, 0, NX*NY*NZ*sizeof(float),
                        (void*)chem1_gpu, 0, NULL, NULL);
    err |= clEnqueueWriteBuffer(cmd_queue, chem2_gpu_mem, CL_TRUE, 0, NX*NY*NZ*sizeof(float),
                        (void*)chem2_gpu, 0, NULL, NULL);
    err |= clEnqueueWriteBuffer(cmd_queue, grid_gpu_mem, CL_TRUE, 0, NX*NY*NZ*sizeof(int),
                        (void*)grid_gpu, 0, NULL, NULL);
    err |= clEnqueueWriteBuffer(cmd_queue, grid1_gpu_mem, CL_TRUE, 0, NX*NY*NZ*sizeof(int),
                        (void*)grid1_gpu, 0, NULL, NULL);
    err |= clEnqueueWriteBuffer(cmd_queue, grid2_gpu_mem, CL_TRUE, 0, NX*NY*NZ*sizeof(int),
                        (void*)grid2_gpu, 0, NULL, NULL);
    err |= clEnqueueWriteBuffer(cmd_queue, grid3_gpu_mem, CL_TRUE, 0, NX*NY*NZ*sizeof(int),
                        (void*)grid3_gpu, 0, NULL, NULL);
    assert(err == CL_SUCCESS);

    {


        size_t global_work_size = cell_no;
        err = clEnqueueNDRangeKernel(cmd_queue, kernel[4], 1, NULL,
                                     &global_work_size, NULL, 0, NULL, NULL);
        assert(err == CL_SUCCESS);
        clFinish(cmd_queue);

    int i_cycle = 0;
    int flag_sort_direction = 0;
    int flag_division = 0;
    int start_ele[MAX_CELL];
        start_ele[0] = 0;
    int iterations = MAX_ITERATIONS;
    int judgedivision = 2;
    for(int iter = 0; iter < iterations; iter++)
	{
         if(iter>=65000 && iter<85000) lx=13.65f-(13.65f-13.00f)/20000*(85000-iter);
         if(iter>=85000 && iter<100000) lx=13.65f+0.1f/15000*(iter-85000);
         if(iter>=100000 && iter<=120000) lx = 13.75f + 1.25f/20000*(iter-100000);
        flag_division = 0;
        for(int i = 0; i < ele_no; i++)
        {
            randx[i] = RNGn(0,1);
            randy[i] = RNGn(0,1);
            randz[i] = RNGn(0,1);
        }
        err = clEnqueueWriteBuffer(cmd_queue, randx_mem, CL_TRUE, 0, max_ele_no*sizeof(float),
                                   (void*)randx, 0, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, randy_mem, CL_TRUE, 0, max_ele_no*sizeof(float),
                                    (void*)randy, 0, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, randz_mem, CL_TRUE, 0, max_ele_no*sizeof(float),
                                    (void*)randz, 0, NULL, NULL);
        err |= clEnqueueWriteBuffer(cmd_queue, ele_per_cell_mem, CL_TRUE, 0, MAX_CELL*sizeof(int),
                                    (void*)ele_per_cell, 0, NULL, NULL);
        assert(err == CL_SUCCESS);
		err = clSetKernelArg(kernel[0],  5, sizeof(int), &cell_no);
		err |= clSetKernelArg(kernel[0], 6, sizeof(int), &ele_no);
        err |= clSetKernelArg(kernel[0], 15, sizeof(float), &lx);
		err |= clSetKernelArg(kernel[0],20, sizeof(int), &i_cycle);
        err |= clSetKernelArg(kernel[0],21, sizeof(cl_mem), &id_is_out_mem);
        err |= clSetKernelArg(kernel[0],29, sizeof(cl_mem), &colevel_mem);
        err |= clSetKernelArg(kernel[0],22, sizeof(int), &iter);
        err |= clSetKernelArg(kernel[0],23, sizeof(int), &flag_division);
		assert(err == CL_SUCCESS);
		global_work_size = ele_no;
        int inner_iter = 1;
        if(i_cycle >= 4) inner_iter = 2;
        for(int i=0; i < inner_iter; i++)
        {
		    err = clEnqueueNDRangeKernel(cmd_queue, kernel[0], 1, NULL, 
									 &global_work_size, NULL, 0, NULL, NULL);
		    assert(err == CL_SUCCESS);
		    clFinish(cmd_queue);
        }

        err = clEnqueueCopyBuffer(cmd_queue, xF_mem, x_mem, 0, 0, max_ele_no*sizeof(float),
                                  0, 0, NULL);
        err |= clEnqueueCopyBuffer(cmd_queue, yF_mem, y_mem, 0, 0, max_ele_no*sizeof(float),
                                  0, 0, NULL);
        err |= clEnqueueCopyBuffer(cmd_queue, zF_mem, z_mem, 0, 0, max_ele_no*sizeof(float),
                                  0, 0, NULL);
		assert(err == CL_SUCCESS);
		clFinish(cmd_queue);

		err = clSetKernelArg(kernel[4],  12, sizeof(int), &ele_no);
        err |= clSetKernelArg(kernel[4],  9, sizeof(float), &lx);
		assert(err == CL_SUCCESS);
        size_t global_work_size = cell_no;
        err = clEnqueueNDRangeKernel(cmd_queue, kernel[4], 1, NULL,
                                     &global_work_size, NULL, 0, NULL, NULL);
        assert(err == CL_SUCCESS);
        clFinish(cmd_queue);

        if(((iter+1)%DIVISION_FR == 0 && judgedivision == 2 && cell_no<32) || (cell_no>=32 && judgedivision == 1 && ((iter+30)*10)%DIVISION_FR == 0 && cell_no < 128))
        {
            err = clEnqueueReadBuffer(cmd_queue, id_mem, CL_TRUE, 0, max_ele_no*sizeof(int),
                              id, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(cmd_queue, x_mem, CL_TRUE, 0, max_ele_no*sizeof(float),
                              x, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(cmd_queue, y_mem, CL_TRUE, 0, max_ele_no*sizeof(float),
                              y, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(cmd_queue, z_mem, CL_TRUE, 0, max_ele_no*sizeof(float),
                              z, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(cmd_queue, ele_type_mem, CL_TRUE, 0, max_ele_no*sizeof(int),
                              ele_type, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(cmd_queue, ele_per_cell_mem, CL_TRUE, 0, cell_no*sizeof(int),
                              ele_per_cell, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(cmd_queue, xc_mem, CL_TRUE, 0, MAX_CELL*sizeof(float),
                              xc, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(cmd_queue, yc_mem, CL_TRUE, 0, MAX_CELL*sizeof(float),
                              yc, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(cmd_queue, zc_mem, CL_TRUE, 0, MAX_CELL*sizeof(float),
                              zc, 0, NULL, NULL);
            assert(err == CL_SUCCESS);

            int tmp_cell_no = cell_no;
            float tmpx0, tmpy0, tmpz0;
            int tmpn0, tmptype0;
            float tmpx[max_ele_no], tmpy[max_ele_no], tmpz[max_ele_no];
            int tmptype[max_ele_no];
            int start, mid;
            int lcell = ele_per_cell[0];
            int tmp_id;
            for(int i = 0; i < cell_no; i++)
            {
                if((i_cycle < 7 && judgedivision==2 && cell_no<32) || (cell_no>=32 && judgedivision==1 && cycle[i]<=0 && ele_per_cell[i]>10))
                {
                lcell = ele_per_cell[i];
                start = start_ele[i];
                mid = start + lcell/2;
                tmpx0 = tmpy0 = tmpz0 = 0.0;
                tmpn0 = tmptype0 = 0;

                tmp_id = id[start];
                for(int j = 0; j < ele_per_cell[i]; j++)
                {
                    int k = start + j;
                    tmpx0 += x[k];
                    tmpy0 += y[k];
                    tmpz0 += z[k];

                    tmpx[j] = x[k];
                    tmpy[j] = y[k];
                    tmpz[j] = z[k];
                    tmptype[j] = ele_type[k];

                    if(id[k] == tmp_id) tmpn0++;
                }

                tmpx0 = tmpx0/tmpn0; //printf("division::cell[%d] tmpx0=%g\n", i, tmpx0);
                tmpy0 = tmpy0/tmpn0;
                tmpz0 = tmpz0/tmpn0;
                    
                int imin;
                float tt;
                int ll = ele_per_cell[i];
                for(int ii=0; ii < ll-1; ii++)
                {
                    imin = ii;
                    if(iter<=260000)
                    {
                        for(int jj=ii+1; jj < ll; jj++)
                        {
                            if(tmpx[imin]>tmpx[jj])
                            {
                                imin = jj;
                            }
                        }
                    }
                    if(imin != ii)
                    {
                        tt = tmpx[imin];
                        tmpx[imin] = tmpx[ii];
                        tmpx[ii] = tt;

                        tt = tmpy[imin];
                        tmpy[imin] = tmpy[ii];
                        tmpy[ii] = tt;

                        tt = tmpz[imin];
                        tmpz[imin] = tmpz[ii];
                        tmpz[ii] = tt;
                    }
                }
                flag_sort_direction++;

                cell_type[tmp_id-1] = cell_type[tmp_cell_no] = 1; 
                id_is_out[tmp_cell_no] = id_is_out[tmp_id-1];
                change_type[tmp_cell_no] = change_type[tmp_id-1];
                change_time[tmp_cell_no] = change_time[tmp_id-1];
                if(iter>60000)
                {
                judgeout[tmp_cell_no] = judgeout[tmp_id-1];
                judgein[tmp_cell_no] = judgein[tmp_id-1];
                }
                if(iter<=60000)
                {
                judgeout[tmp_cell_no] = 50;
                judgein[tmp_cell_no] = 50;
                judgeout[tmp_id-1] = 50;
                judgein[tmp_id-1] = 50;
                }
                colevel[tmp_cell_no] = colevel[tmp_id-1];
                mother_id[tmp_id-1]=id_is_out[tmp_id-1];
                mother_id[tmp_cell_no]=id_is_out[tmp_id-1];
                if(ele_per_cell[i] == 40)
                {
                    cycle[tmp_id-1]=10+RNG()*20;
                    cycle[tmp_cell_no]=10+RNG()*20;
                }
                if(ele_per_cell[i] == 20)
                {
                    cycle[tmp_id-1]=10000;
                    cycle[tmp_cell_no]=10000;
                }
                signal[tmp_cell_no][0] = signal[tmp_id-1][0];
                signal[tmp_cell_no][1] = signal[tmp_id-1][1];
                signal[tmp_cell_no][2] = signal[tmp_id-1][2];
                signal[tmp_cell_no][3] = signal[tmp_id-1][3];
                signal[tmp_cell_no][4] = signal[tmp_id-1][4];
                signal[tmp_cell_no][5] = signal[tmp_id-1][5];
                signal[tmp_cell_no][6] = signal[tmp_id-1][6];

                int tmp_ele_per_cell = ele_per_cell[i];
                ele_per_cell[i] = 0;
                ele_per_cell[tmp_cell_no] = 0;
                int tmp_count1 = 0, tmp_count2 = 0;
                for(int j = 0; j < tmp_ele_per_cell; j++)
                {
                    if(tmpx[j] > tmpx0) tmp_count1++;
                    else tmp_count2++;
                }
                int tt_half = tmp_ele_per_cell/2;
                start_ele[tmp_cell_no] = mid;
                for(int j = 0; j < tt_half; j++)
                {
                    x[start] = tmpx[j];
                    y[start] = tmpy[j];
                    z[start] = tmpz[j];
                    id[start] = i+1;
                    ele_type[start] = tmptype[j];
                    ele_per_cell[i]++;
                    start++;

                    x[mid] = tmpx[tt_half+j];
                    y[mid] = tmpy[tt_half+j];
                    z[mid] = tmpz[tt_half+j];
                    id[mid] = tmp_cell_no+1;
                    ele_type[mid] = tmptype[tt_half+j];
                    ele_per_cell[tmp_cell_no]++;
                    mid++;
                }

                tmp_cell_no++;
                }
            else if(cell_no>=32 && judgedivision==1 && cycle[i]>0)
            {
            cycle[i]=cycle[i]-1;
            }
            }
            cell_no = tmp_cell_no;

            err = clEnqueueWriteBuffer(cmd_queue, id_mem, CL_TRUE, 0, max_ele_no*sizeof(int),
                            (void*)id, 0, NULL, NULL);
            err |= clEnqueueWriteBuffer(cmd_queue, x_mem, CL_TRUE, 0, max_ele_no*sizeof(float),
                            (void*)x, 0, NULL, NULL);
            err |= clEnqueueWriteBuffer(cmd_queue, y_mem, CL_TRUE, 0, max_ele_no*sizeof(float),
                            (void*)y, 0, NULL, NULL);
            err |= clEnqueueWriteBuffer(cmd_queue, z_mem, CL_TRUE, 0, max_ele_no*sizeof(float),
                            (void*)z, 0, NULL, NULL);
            err |= clEnqueueWriteBuffer(cmd_queue, cell_type_mem, CL_TRUE, 0, cell_no*sizeof(int),
                            (void*)cell_type, 0, NULL, NULL);
            err |= clEnqueueWriteBuffer(cmd_queue, ele_type_mem, CL_TRUE, 0, max_ele_no*sizeof(int),
                            (void*)ele_type, 0, NULL, NULL);
            err |= clEnqueueWriteBuffer(cmd_queue, ele_per_cell_mem, CL_TRUE, 0, MAX_CELL*sizeof(int),
                            (void*)ele_per_cell, 0, NULL, NULL);
            assert(err == CL_SUCCESS);

        if(i_cycle<6) i_cycle++;
        if(i_cycle == 5) i_cycle=6;


        flag_division = 1;
        err = clSetKernelArg(kernel[0],  5, sizeof(int), &cell_no);
        err |= clSetKernelArg(kernel[0], 6, sizeof(int), &ele_no);
        err |= clSetKernelArg(kernel[0],20, sizeof(int), &i_cycle);
        err |= clSetKernelArg(kernel[0],21, sizeof(cl_mem), &id_is_out_mem);
        err |= clSetKernelArg(kernel[0],29, sizeof(cl_mem), &colevel_mem);
        err |= clSetKernelArg(kernel[0],22, sizeof(int), &iter);
        err |= clSetKernelArg(kernel[0],23, sizeof(int), &flag_division);
        assert(err == CL_SUCCESS);
        global_work_size = ele_no;
        int inner_iter = 10000;
        for(int i=0; i < inner_iter; i++)
        {
            err = clEnqueueNDRangeKernel(cmd_queue, kernel[0], 1, NULL,
                                     &global_work_size, NULL, 0, NULL, NULL);
            assert(err == CL_SUCCESS);
            clFinish(cmd_queue);
        }

        err = clEnqueueCopyBuffer(cmd_queue, xF_mem, x_mem, 0, 0, max_ele_no*sizeof(float),
                                  0, 0, NULL);
        err |= clEnqueueCopyBuffer(cmd_queue, yF_mem, y_mem, 0, 0, max_ele_no*sizeof(float),
                                  0, 0, NULL);
        err |= clEnqueueCopyBuffer(cmd_queue, zF_mem, z_mem, 0, 0, max_ele_no*sizeof(float),
                                  0, 0, NULL);
        assert(err == CL_SUCCESS);
        clFinish(cmd_queue);


        }

        
        err = clSetKernelArg(kernel[2],  7, sizeof(int), &cell_no);
        err |= clSetKernelArg(kernel[2], 0, sizeof(cl_mem), &xc_mem);
        err |= clSetKernelArg(kernel[2], 3, sizeof(float), &lx);
        err |= clSetKernelArg(kernel[2], 1, sizeof(cl_mem), &yc_mem);
        err |= clSetKernelArg(kernel[2], 2, sizeof(cl_mem), &zc_mem);
        assert(err == CL_SUCCESS);
        err = clEnqueueWriteBuffer(cmd_queue, id_is_out_mem, CL_TRUE, 0, cell_no*sizeof(int),
                                   (void*)id_is_out, 0, NULL, NULL);
        assert(err == CL_SUCCESS);
        err = clEnqueueWriteBuffer(cmd_queue, colevel_mem, CL_TRUE, 0, cell_no*sizeof(int),
                                   (void*)colevel, 0, NULL, NULL);
        assert(err == CL_SUCCESS);
        global_work_size = cell_no;
        err = clEnqueueNDRangeKernel(cmd_queue, kernel[2], 1, NULL,
                                     &global_work_size, NULL, 0, NULL, NULL);
        assert(err == CL_SUCCESS);
        clFinish(cmd_queue);
        
        err = clEnqueueReadBuffer(cmd_queue, flag_contact_mem, CL_TRUE, 0, cell_no*sizeof(int),
                                  flag_contact, 0, NULL, NULL);
        assert(err == CL_SUCCESS);


        if(i_cycle >= 0 && iter%10 == 0 ) // decrease frequency
        {
            err = clEnqueueReadBuffer(cmd_queue, flag_outer_mem, CL_TRUE, 0, MAX_CELL*sizeof(int),
                                  flag_outer, 0, NULL, NULL);
            assert(err == CL_SUCCESS);
            if(iter <= 260000){
            int ii0 = 1, jj0 = 4;
            for(int i=0; i<ii0; i++)
            {
                for(int j=0; j<jj0; j++){
                signal_update(cell_no, dt, iter+1, flag_outer, ele_per_cell, i_cycle, flag_contact, judgeout, judgein, colevel, mother_id);
                }
                int * countcon = (int *)malloc(MAX_CELL*sizeof(int));
                float * signal_noise = (float *)malloc(MAX_CELL*sizeof(float));
                float * dis_save = (float *)malloc(MAX_CELL*sizeof(float));
                int * contactnum = (int *)malloc(MAX_CELL*sizeof(int));
                int * contactin = (int *)malloc(MAX_CELL*sizeof(int));
                int * contactout = (int *)malloc(MAX_CELL*sizeof(int));
                int tmpid;
                int tmpidio;
                float sum_dis;
                int cid;
                float contis;
                float judgedis;
                int receive_noise;
                for(int i=0; i<cell_no; i++){
                    signal[i][7] = 0;
                    signal_noise[i] = RNG();
                }
                
                for(int i=0; i<cell_no; i++){
                    contactin[i] = 0;
                    contactout[i] = 0;
                    for(int j=0; j<cell_no; j++){
                        if(ele_per_cell[i]+ele_per_cell[j]==20) judgedis=1.0*1.7*3.79f;
                        if(ele_per_cell[i]+ele_per_cell[j]==30) judgedis=1.0*1.7*4.189f;
                        if(ele_per_cell[i]+ele_per_cell[j]==40) judgedis=1.0*1.7*4.588f;
                        if(ele_per_cell[i]+ele_per_cell[j]==50) judgedis=1.0*1.7*4.9879f;
                        if(ele_per_cell[i]+ele_per_cell[j]==60) judgedis=1.0*1.7*5.387f;
                        if(ele_per_cell[i]+ele_per_cell[j]==80) judgedis=1.0*1.7*6.1859f;
                        if(ele_per_cell[i]+ele_per_cell[j]==160) judgedis=1.0*1.7*8.93f;
                        if(ele_per_cell[i]+ele_per_cell[j]>=320) judgedis=1.0*1.7*10.83f;
                        contis = distance3(xc[i], yc[i], zc[i], xc[j], yc[j], zc[j]);
                        if(contis < 1.0*judgedis && contis > -0.0001f){
                            if(id_is_out[j] == 2) contactout[i]++;
                            else contactin[i]++;
                        }
                    }
                }
                
                for(int i=0; i<cell_no; i++){
                    if(id_is_out[i] != 2 || iter<50000){
                        tmpidio = 0;
                        tmpid = 0;
                        sum_dis =0.0f;
                        for(int j=0; j<cell_no; j++){
                            if(ele_per_cell[i]+ele_per_cell[j]==20) judgedis=1.5*1.7*3.79f;
                            if(ele_per_cell[i]+ele_per_cell[j]==30) judgedis=1.5*1.7*4.189f;
                            if(ele_per_cell[i]+ele_per_cell[j]==40) judgedis=1.5*1.7*4.588f;
                            if(ele_per_cell[i]+ele_per_cell[j]==50) judgedis=1.5*1.7*4.9879f;
                            if(ele_per_cell[i]+ele_per_cell[j]==60) judgedis=0.5*1.7*5.387f;
                            if(ele_per_cell[i]+ele_per_cell[j]==80) judgedis=0.5*1.7*6.1859f;
                            if(ele_per_cell[i]+ele_per_cell[j]==160) judgedis=0.5*1.7*8.93f;
                            if(ele_per_cell[i]+ele_per_cell[j]>=320) judgedis=0.5*1.7*10.83f;
                            contis = distance3(xc[i], yc[i], zc[i], xc[j], yc[j], zc[j]);
                            if(contis < 1.0*judgedis && contis > -0.0001f){
                                tmpidio++;
                                if(id_is_out[j] != 2 || iter < 50000){
                                    if(iter > 50000) receive_noise = 1.0 + 0.0*RNG();
                                    if(iter <= 50000) receive_noise = 1.0;
                                    countcon[tmpid] = j;
                                    dis_save[tmpid] = receive_noise*1.0f;
                                    if(iter>=120000){
                                        if(contactout[j] > 2 && contactin[j] < 9) dis_save[tmpid] = receive_noise*1.8f;
                                        if(contactout[j] < 2 && contactin[j] < 9) dis_save[tmpid] = receive_noise*1.4f;
                                    }
                                    if(iter>=95000 && iter<120000){
                                        if(contactin[j] < 6) dis_save[tmpid] = receive_noise*1.4f;
                                    }
                                    if(iter>=50000 && iter<95000){
                                        if(contactin[j] < 5) dis_save[tmpid] = receive_noise*1.4f;
                                    }
                                    sum_dis += dis_save[tmpid];
                                    tmpid++;
                                    
                                }
                            }
                        }
                        if(iter == 150000 || iter == 65000 || iter == 75000 || iter == 95000){
                        printf("N_with_TE:%d, N_without_TE:%d \n",contactout[i], contactin[i]);
                        }
                        contactnum[i] = tmpid;
                        if(tmpid > 0){
                            for(int k=0; k<tmpid; k++){
                                cid = countcon[k];
                                signal[cid][7] = signal[cid][7] + (0.9 + 0.2*signal_noise[cid])*dis_save[k]*signal[i][6]/sum_dis;
                            }
                        }
                    }
                }
                
            }
                if(iter >= 65000){
                    float xcc = 0.0f, ycc = 0.0f, zcc = 0.0f, xyzcc;
                    int idin = 0;
                    for(int i=0; i<cell_no; i++){
                        if(id_is_out[i] != 2){
                            xcc = xcc+xc[i];
                            ycc = ycc+yc[i];
                            zcc = zcc+zc[i];
                            idin++;
                        }
                    }
                    xcc = xcc/idin;
                    ycc = ycc/idin;
                    zcc = zcc/idin;
                    xyzcc = distance3(xcc, ycc, zcc, 0.0f, 0.0f, 0.0f);
                    xcc = -xcc/xyzcc;
                    ycc = -ycc/xyzcc;
                    zcc = -zcc/xyzcc;
                    err = clSetKernelArg(kernel[0],  30, sizeof(float), &xcc);
                    err |= clSetKernelArg(kernel[0], 31, sizeof(float), &ycc);
                    err |= clSetKernelArg(kernel[0], 32, sizeof(float), &zcc);
                    assert(err == CL_SUCCESS);
                }

            err = clEnqueueWriteBuffer(cmd_queue, id_is_out_mem, CL_TRUE, 0, cell_no*sizeof(int),
                            (void*)id_is_out, 0, NULL, NULL);
            assert(err == CL_SUCCESS);
            err = clEnqueueWriteBuffer(cmd_queue, colevel_mem, CL_TRUE, 0, cell_no*sizeof(int),
                            (void*)colevel, 0, NULL, NULL);
            assert(err == CL_SUCCESS);


            if(flag_print == 1)
            {
                for(int i=0; i<cell_no; i++)
                {
                    printf("cell[%d]:outer ele %d\n",i,flag_outer[i]);
                }
            }
            }
            if(iter == 75000){
                for(int i=0; i<cell_no; i++){
                    cycle[i]=RNG()*20;
                }
                judgedivision=1;
            }
        }

        // output
        if(iter%OUTPUT_FR == 0)
        {
		    err = clEnqueueReadBuffer(cmd_queue, x_mem, CL_TRUE, 0, max_ele_no*sizeof(float), 
								  x, 0, NULL, NULL);
		    err |= clEnqueueReadBuffer(cmd_queue, y_mem, CL_TRUE, 0, max_ele_no*sizeof(float), 
								  y, 0, NULL, NULL);
		    err |= clEnqueueReadBuffer(cmd_queue, z_mem, CL_TRUE, 0, max_ele_no*sizeof(float), 
								  z, 0, NULL, NULL);
		    err |= clEnqueueReadBuffer(cmd_queue, id_mem, CL_TRUE, 0, max_ele_no*sizeof(int), 
								  id, 0, NULL, NULL);
		    err |= clEnqueueReadBuffer(cmd_queue, ele_type_mem, CL_TRUE, 0, max_ele_no*sizeof(int), 
								  ele_type, 0, NULL, NULL);
		    err |= clEnqueueReadBuffer(cmd_queue, ele_per_cell_mem, CL_TRUE, 0, MAX_CELL*sizeof(int), 
								  ele_per_cell, 0, NULL, NULL);
		    err |= clEnqueueReadBuffer(cmd_queue, chem1_gpu_mem, CL_TRUE, 0, NX*NY*NZ*sizeof(float), 
								  chem1_gpu, 0, NULL, NULL);
		    err |= clEnqueueReadBuffer(cmd_queue, chem2_gpu_mem, CL_TRUE, 0, NX*NY*NZ*sizeof(float), 
								  chem2_gpu, 0, NULL, NULL);
		    assert(err == CL_SUCCESS);
            err = clEnqueueReadBuffer(cmd_queue, xc_mem, CL_TRUE, 0, cell_no*sizeof(float),
                                  xc, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(cmd_queue, yc_mem, CL_TRUE, 0, cell_no*sizeof(float),
                                  yc, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(cmd_queue, zc_mem, CL_TRUE, 0, cell_no*sizeof(float),
                                  zc, 0, NULL, NULL);
            err |= clEnqueueReadBuffer(cmd_queue, ele_per_cell_mem, CL_TRUE, 0, MAX_CELL*sizeof(int),
                                  ele_per_cell, 0, NULL, NULL);
            assert(err == CL_SUCCESS);

            output(x, y, z, id, cell_type, ele_type, cell_no, ele_no, ele_per_cell, iter, chem1_gpu, chem2_gpu, xc, yc, zc, flag_outer, id_is_out, signal, colevel, id_is_f);

            if(iter >= DIVISION_FR*26)
            {
                printf("END::iter=%d i_cycle=%d cell_no=%d ele_no=%d\n", iter, i_cycle, cell_no, ele_no);
                exit(0);
            }
        }



        // statistics
        if(iter%1000 == 0)
        {
            int count0 = 0, count1 = 0, count2 = 0, count3 = 0;
            int c_position_out = 0, c_position_in = 0, c_in_in = 0, c_in_in1 = 0, c_in_in2 = 0, c_out_out = 0;
            for(int i = 0; i < cell_no; i++)
            {
                count0++;
                if(id_is_out[i] != 2) count1++;
                if(id_is_out[i] == 0) count2++;
                if(id_is_out[i] == 2) count3++;
                if(flag_outer[i] >= 2) c_position_out++;
                if(flag_outer[i] == 0) c_position_in++;
                if(flag_outer[i] >= 2 && id_is_out[i] == 2) c_out_out++;
                if(flag_outer[i] == 0 && id_is_out[i] == 0) c_in_in++;
                if(flag_outer[i] == 0 && id_is_out[i] == 1) c_in_in1++;
                if(flag_outer[i] == 0 && id_is_out[i] == 4) c_in_in2++;
            }
            printf("%d %6d %4d   %4d %4d %4d %4d  %4d %4d\n", i_cycle, iter, count0, c_position_in, c_in_in, c_in_in1, c_in_in2, c_position_out, c_out_out);
        }

    }
       
		
	}
	
#pragma mark Teardown

	{
		clReleaseMemObject(x_mem);
		clReleaseMemObject(y_mem);
		clReleaseMemObject(z_mem);
		clReleaseMemObject(xc_mem);
		clReleaseMemObject(yc_mem);
		clReleaseMemObject(zc_mem);
		clReleaseMemObject(id_mem);
		clReleaseMemObject(cell_type_mem);
		clReleaseMemObject(ele_type_mem);
		clReleaseMemObject(xF_mem);
		clReleaseMemObject(yF_mem);
		clReleaseMemObject(zF_mem);
		clReleaseMemObject(ele_per_cell_mem);
		clReleaseMemObject(flag_mem);
		clReleaseMemObject(chem1_gpu_mem);
		clReleaseMemObject(chem2_gpu_mem);
		clReleaseMemObject(chem11_gpu_mem);
		clReleaseMemObject(chem22_gpu_mem);
		clReleaseMemObject(grid_gpu_mem);
		clReleaseMemObject(grid1_gpu_mem);
		clReleaseMemObject(grid2_gpu_mem);
		clReleaseMemObject(grid3_gpu_mem);
		
		clReleaseCommandQueue(cmd_queue);
		clReleaseContext(context);
	}
	return CL_SUCCESS;
}


void  cell_in_block(float * x, float * y, float * z, int * id, int * cell_type, int * ele_type, int cell_no, int * ele_per_cell, int * judgeout, int * judgein)
{
    int i, j, k;
    int ni, nj, nk;

    for(i = 0; i < NX; i++)
        for(j = 0; j < NY; j++)
            for(k = 0; k < NZ; k++)
    {
        grid[i][j][k] = 0;
        grid1[i][j][k] = 0;
        grid2[i][j][k] = 0;
        grid3[i][j][k] = 0;
    }

    for(i = 0; i < cell_no; i++)
    {
        for(j = 0; j < ele_per_cell[i]; j++)
        {
            k = i*MAX_ELE + j;
            ni = (int)(1.0*x[k]/(1.0*DX)); if(ni < 0) ni = 0; if(ni >= NX) ni = NX-1;
            nj = (int)(1.0*y[k]/(1.0*DY)); if(nj < 0) nj = 0; if(nj >= NY) nj = NY-1;
            nk = (int)(1.0*z[k]/(1.0*DZ)); if(nk < 0) nk = 0; if(nk >= NZ) nk = NZ-1;

            grid[ni][nj][nk]++;
            if(cell_type[i] == 1)
            {
                grid1[ni][nj][nk]++;
            }
            else if(cell_type[i] == 2)
            {
                grid2[ni][nj][nk]++;
            }
            else if(cell_type[i] == 3)
            {
                grid3[ni][nj][nk]++;
            }
            else
            {
                printf("ERROR::cell_in_block(), wrong cell type = %d\n", cell_type[i]);
                exit(0);
            }
        }
    } 

}

void  initial_chem_paras(float * x, float * y, float * z, int * id, int * cell_type, int * ele_type, int cell_no, int * ele_per_cell)
{
    int i, j, k;

    for(i = 0; i < NX; i++)
        for(j = 0; j < NY; j++)
            for(k = 0; k < NZ; k++)
    {
        grid[i][j][k] = 0;
        chem1[i][j][k] = 0.0;
        chem2[i][j][k] = 0.0;
    }
    return;
}


void signal_init()
{
    for(int i=0; i < MAX_CELL; i++)
    {
        signal[i][0] = 0.0f;
        signal[i][0] = 0.0f;
        signal[i][0] = 1.0f;
        signal[i][0] = 1.0f;
        signal[i][2] = 1.0f;
        signal[i][2] = 1.0f;
        signal[i][3] = 1.0f;
        signal[i][3] = 1.0f;
    }
}

void signal_update(int cell_no, float dt, int iter, int * flag_outer, int * ele_per_cell, int i_cycle, int *flag_contact, int * judgeout, int * judgein, int * colevel, int * mother_id)
{
    dt = 4*dt;
    float dy[10];
    
    float ax = 1.0f;
    float ay = 1.0f;
    int n = 4;
    float theta = 0.5f;
    float delta = 1.0f;
    
    float Ix = 0.4f;
    float Iy = 0.4f;
    
    float k = 0.8f; // 1.7; 1.4;
    float lambda1 = 0.4f;// 0.4;
    
    float bx = 2.0f;
    float by = 0.7f;
    float bg = 0.88f;
    float ratiooo=0.03f;
    
    float vsg1 = 1.2f;
    float vsg2 = 1.0f;
    float vsn1 = 0.856f;
    float vsn2 = 1.0f;
    float vsfr1 = 2.8f;
    float vsfr2 = 2.8f;
    float vex = 0.0f;
    float vsf = 0.6f;
    float va = 20.0f;
    float vi = 3.3f;
    float kdg = 1.0f;
    float kdn = 1.0f;
    float kdfr = 1.0f;
    float kdf = 0.077f;
    float kag1 = 0.28f;
    float kag2 = 0.55f;
    float kan = 0.55f;
    float kafr = 0.5f;
    float kaf = 5.0f;
    float kig = 0.8*2.0f;
    float kin1 = 0.28f;
    float kin2 = 0.802*2.0f;
    float kifr = 0.5f;
    float ka = 0.7f;
    float ki = 0.7f;
    float kd = 2.0f;
    kin1 = 0.28f;
    kag1 = 0.28f;
    int pr = 3;
    int ps = 4;
    int pq = 4;
    int pu = 3;
    int pv = 4;
    int pw = 4;
    int px = 1;
    int py = 1;
    int pz = 4;
    
    float theta_n = pow(theta, n);
    float kag1_r = pow(kag1, pr);
    float kag2_s = pow(kag2, ps);
    float kig_q = pow(kig, pq);
    float kin1_u = pow(kin1, pu);
    float kan_v = pow(kan, pv);
    float kin2_w = pow(kin2, pw);
    float kifr_x = pow(kifr, px);
    float kafr_y = pow(kafr, py);
    float kaf_z = pow(kaf, pz);
    
    float s0n, s1n, s2q, s2v, s2x, s2z, s3s, s3w, s3y, s5r, s5u, flag1, flag2, tmp, tmp1;
    
    float tmp_diff = 0.6f;
    
    for(int i=0; i<cell_no; i++)
    {
        s0n = pow(signal[i][0],n);
        s1n = pow(signal[i][1],n);
        s2q = pow(signal[i][2],pq);
        s2v = pow(signal[i][2],pv);
        s2x = pow(signal[i][2],px);
        s2z = pow(signal[i][2],pz);
        s3s = pow(signal[i][3],ps);
        s3w = pow(signal[i][3],pw);
        s3y = pow(signal[i][3],py);
        s5r = pow(signal[i][5],pr);
        s5u = pow(signal[i][5],pu);
        float fr = signal[i][4];
        float erk = signal[i][5];
        float fp = signal[i][7];
        tmp = 1.0f*flag_outer[i]/ele_per_cell[i];
        tmp = 1.0f - tmp;
        tmp1 = flag_contact[i];
        if(tmp < 0.0f) tmp = 0.0f;
        bx = 2.0f - 1.5f*tmp;
        tmp1 = (tmp1-5.0f)/5.0f;
        if(tmp1>1.0f) tmp1=1.0f;
        if(tmp1<0.0f) tmp1=0.0f;
        tmp1=1.0f-tmp1;
        float ri;
        if(iter>=280000)
        {
           bg = bg - 0.0*tmp1;
        }
        if(iter>160000)
        {
            vsn1=vsn1*1.0f;
            vsg1=vsg1*1.0f;
        }
        
        if(iter<=265000)
        {
            dy[0] = k*(bx + ax*s0n/(theta_n + s0n))*(Ix + (1-Ix)*theta_n/(theta_n + s1n)) - signal[i][0];
            dy[1] = k*(by + ay*s1n/(theta_n + s1n))*(Iy + (1-Iy)*theta_n/(theta_n + s0n)) - delta*signal[i][1];
        }
        
        vsn1 = 0.856f;
        vsg1 = 0.856f;
        vsn2 = 1.0f;
        vsg2 = 1.0f;
        kan = 0.55f;
        kag2 = 0.55f;
        kin2 = 0.802*2.0f;
        kig = 0.8*2.0f;
        ps = 4;
        pv = 4;
        pw = 4;
        pq = 4;
        kdn = 1.0f;
        kdg = 1.0f;
        
        if((id_is_out[i]!=2 || iter<=60000) && iter >= 0)
        {
            ri = 1.0f;
            dy[2] = 1.15*0.040752 - 15*0.2*0.002*vsn1*ri*s5u*kin1_u/(kin1_u + ri*s5u) + 0.321*vsn2*s2v/(kan_v + s2v)*kin2_w/(kin2_w + s3w) - 0.2*kdn*signal[i][2];
            dy[3] = 0.04103217 + 15*0.2*0.0112*vsg1*ri*s5r/(kag1_r + ri*s5r) + 0.321*vsg2*s3s/(kag2_s + s3s)*kig_q/(kig_q + s2q) - 0.2*kdg*signal[i][3];
            dy[4] = vsfr1*kifr_x/(kifr_x + s2x) + vsfr2*s3y/(kafr_y + s3y) - kdfr*signal[i][4];
            dy[5] = va*fr*fp/(kd + fp)*(1 - erk)/(ka + 1 - erk) - vi*erk/(ki + erk);
            dy[6] = 1.2*vex + vsf*s2z/(kaf_z + s2z) - kdf*signal[i][6];
        }
        
        if(id_is_out[i]!=2 && iter>60000)
        {
            ri = 1.0 - 1.0f*(iter-60000.0)/70000.0;
            ri = 1.0f;
            if(iter < FGF_on_time) ri = 0.0f;
            if(iter >= FGF_off_time) ri = 0.0f;
            dy[2] = 1.15*0.040752 - 20*0.2*0.002*vsn1*ri*s5u*kin1_u/(kin1_u + ri*s5u) + 0.321*vsn2*s2v/(kan_v + s2v)*kin2_w/(kin2_w + s3w) - 0.2*kdn*signal[i][2];
            dy[3] = 0.04103217 + 20*0.2*0.0112*vsg1*ri*s5r/(kag1_r + ri*s5r) + 0.321*vsg2*s3s/(kag2_s + s3s)*kig_q/(kig_q + s2q) - 0.2*kdg*signal[i][3];
            dy[4] = vsfr1*kifr_x/(kifr_x + s2x) + vsfr2*s3y/(kafr_y + s3y) - kdfr*signal[i][4];
            dy[5] = va*fr*fp/(kd + fp)*(1 - erk)/(ka + 1 - erk) - vi*erk/(ki + erk);
            dy[6] = 1.2*vex + vsf*s2z/(kaf_z + s2z) - kdf*signal[i][6];
        }
        
        if(iter<=265000)
        signal[i][0] = signal[i][0] +  dy[0]*dt  + lambda1*sqrt(dt)*RNGn(0,1)*signal[i][0];
        signal[i][0] = (signal[i][0] > 0.0f)? signal[i][0] : 0.0f;
        signal[i][1] = signal[i][1] +  dy[1]*dt  + lambda1*sqrt(dt)*RNGn(0,1)*signal[i][1];
        signal[i][1] = (signal[i][1] > 0.0f)? signal[i][1] : 0.0f;
        if(iter<260000)
        {
            float noise_level;
            if(iter < 50000) noise_level = 0.0f;
            if(iter >= 50000) noise_level = 1.0f*(iter-50000.0f)/60000.0f;
            if(noise_level < 0.0f) noise_level = 0.0f;
            if(noise_level > 1.0f) noise_level = 1.0f;
            signal[i][2] = signal[i][2] +  ratiooo*dy[2]*dt  + noise_level*0.018*sqrt(ratiooo*dt)*RNGn(0,1)*signal[i][2];
            signal[i][2] = (signal[i][2] > 0.0f)? signal[i][2] : 0.0f;
            signal[i][3] = signal[i][3] +  ratiooo*dy[3]*dt  + noise_level*0.018*sqrt(ratiooo*dt)*RNGn(0,1)*signal[i][3];
            signal[i][3] = (signal[i][3] > 0.0f)? signal[i][3] : 0.0f;
            signal[i][4] = signal[i][4] +  ratiooo*dy[4]*dt;
            signal[i][4] = (signal[i][4] > 0.0f)? signal[i][4] : 0.0f;
            signal[i][5] = signal[i][5] +  ratiooo*dy[5]*dt;
            signal[i][5] = (signal[i][5] > 0.0f)? signal[i][5] : 0.0f;
            signal[i][6] = signal[i][6] +  ratiooo*dy[6]*dt;
            signal[i][6] = (signal[i][6] > 0.0f)? signal[i][6] : 0.0f;
        }
        flag1 = max(0.0f, min(tmp_diff, signal[i][0] - signal[i][1]));
        flag2 = max(0.0f, min(tmp_diff, signal[i][1] - signal[i][0]));
        if(iter>0)
        {
        if(tmp<0.9f && i_cycle>=3)
        {
            if(judgeout[i]<50) judgeout[i] = judgeout[i]+1;
            if(judgein[i]>0) judgein[i] = judgein[i]-1;
        }
        if(tmp>=0.9f && i_cycle>=3)
        {
            if(judgeout[i]>0) judgeout[i] = judgeout[i]-1;
            if(judgein[i]<50) judgein[i] = judgein[i]+1;
        }
        if(id_is_out[i] != 2 && judgeout[i] == 50 && iter<51000)    id_is_out[i]=2;
        if(id_is_out[i] == 2 && judgein[i] == 50 && iter<51000)    id_is_out[i]=1;
        if(i_cycle>=3 && judgeout[i]==50) judgein[i]=0;
        if(i_cycle>=3 && judgein[i]==50) judgeout[i]=0;
            if(flag1 == tmp_diff)
            {
                colevel[i] = 2; // out
            }
            else if(flag2 == tmp_diff)
            {
                colevel[i] = 0;// in
            }
            else
            {
                colevel[i] = 1; // undecided
            }

        }
        
        float k1 = 2.624f;
        float b1 = -0.1937f;
        float k2 = 0.5248f;
        float b2 = 0.08836f;
        if(iter>0 && id_is_out[i] != 2 && iter>50000)
        {
            if(signal[i][2]>signal[i][3]*k1+b1 && signal[i][2]>signal[i][3]*k2+b2) id_is_out[i]=0;
            if(signal[i][2]<signal[i][3]*k1+b1 && signal[i][2]<signal[i][3]*k2+b2) id_is_out[i]=4;
            if(signal[i][2]>signal[i][3]*k1+b1 && signal[i][2]<signal[i][3]*k2+b2) id_is_out[i]=3;
            if(signal[i][2]<signal[i][3]*k1+b1 && signal[i][2]>signal[i][3]*k2+b2) id_is_out[i]=1;
        }
        if(signal[i][4]<1.5 && signal[i][6]<0.03) id_is_f[i]=0;
        if(signal[i][4]<1.5 && signal[i][6]>0.03) id_is_f[i]=1;
        if(signal[i][4]>1.5 && signal[i][6]<0.03) id_is_f[i]=2;
        if(signal[i][4]>1.5 && signal[i][6]>0.03) id_is_f[i]=3;
        
    }
}

int main (int argc, const char * argv[]) {

    int i, j;
   
    char InName[] = "Input/IC.txt";
    int max_ele_no = MAX_CELL*MAX_ELE; //DEFAULT_MAX_ELE_NO;
    int ele_no = 0;
    int cell_no = 0; 
    float dt = 0.1;

	// Allocate some memory and a place for the results
	float * x = (float *)malloc(max_ele_no*sizeof(float));
	float * y = (float *)malloc(max_ele_no*sizeof(float));
	float * z = (float *)malloc(max_ele_no*sizeof(float));
	int * id = (int *)malloc(max_ele_no*sizeof(int));
	int * ele_type = (int *)malloc(max_ele_no*sizeof(int));
	int * cell_type = (int *)malloc(MAX_CELL*sizeof(int));
	int * ele_per_cell = (int *)malloc(MAX_CELL*sizeof(int));

    chem1 = (float ***)malloc(NX*sizeof(float**));
    for(i = 0; i < NX; i++)
    {
        chem1[i] = (float **)malloc(NY*sizeof(float*));
        for(j = 0; j < NY; j++)
            chem1[i][j] = (float *)malloc(NZ*sizeof(float));
    }
    chem2 = (float ***)malloc(NX*sizeof(float**));
    for(i = 0; i < NX; i++)
    {
        chem2[i] = (float **)malloc(NY*sizeof(float*));
        for(j = 0; j < NY; j++)
            chem2[i][j] = (float *)malloc(NZ*sizeof(float));
    }
    chem_diff = (float ***)malloc(NX*sizeof(float**));
    for(i = 0; i < NX; i++)
    {
        chem_diff[i] = (float **)malloc(NY*sizeof(float*));
        for(j = 0; j < NY; j++)
            chem_diff[i][j] = (float *)malloc(NZ*sizeof(float));
    }
    grid = (int ***)malloc(NX*sizeof(int**));
    for(i = 0; i < NX; i++)
    {
        grid[i] = (int **)malloc(NY*sizeof(int*));
        for(j = 0; j < NY; j++)
            grid[i][j] = (int *)malloc(NZ*sizeof(int));
    }
    grid1 = (int ***)malloc(NX*sizeof(int**));
    for(i = 0; i < NX; i++)
    {
        grid1[i] = (int **)malloc(NY*sizeof(int*));
        for(j = 0; j < NY; j++)
            grid1[i][j] = (int *)malloc(NZ*sizeof(int));
    }
    grid2 = (int ***)malloc(NX*sizeof(int**));
    for(i = 0; i < NX; i++)
    {
        grid2[i] = (int **)malloc(NY*sizeof(int*));
        for(j = 0; j < NY; j++)
            grid2[i][j] = (int *)malloc(NZ*sizeof(int));
    }
    grid3 = (int ***)malloc(NX*sizeof(int**));
    for(i = 0; i < NX; i++)
    {
        grid3[i] = (int **)malloc(NY*sizeof(int*));
        for(j = 0; j < NY; j++)
            grid3[i][j] = (int *)malloc(NZ*sizeof(int));
    }
  
    signal = (float **)malloc(MAX_CELL*sizeof(float*));
    for(i = 0; i < MAX_CELL; i++)
    {
        signal[i] = (float *)malloc(10*sizeof(float));
    }
    id_is_out = (int *)malloc(MAX_CELL*sizeof(int));
    id_is_f = (int *)malloc(MAX_CELL*sizeof(int));
    for(i=0;i<MAX_CELL;i++)
    {
        id_is_out[i] = 2; // out
        id_is_f[i] = 0;
        for(j=0;j<10;j++)
        {
            signal[i][j] = 0.0f;
//            if(j==2 || j==3) signal[i][j] =  0+0.01*RNG();
            if(j==2 || j==3) signal[i][j] =  0.5;
            if(j==4) signal[i][j] = 2.8;
            if(j==5) signal[i][j] = 0.25;
            if(j==6) signal[i][j] = 0.07;
            
        }
    }

    InitialReader(InName, x, y, z, id, cell_type, ele_type, &cell_no, &ele_no, ele_per_cell);

    // INIT::Fill in the values
    float tmp_r, theta, phi;
    for(i=0;i<1280;i++)
    {
        tmp_r = 10.0f*RNG();
		theta = PI*RNG();
        phi = 2.0f*PI*RNG();
        x[i] = tmp_r*sin(theta)*cos(phi);
        y[i] = tmp_r*sin(theta)*sin(phi);
        z[i] = tmp_r*cos(theta);
        id[i] = 1;
        ele_type[i] = 1;
        cell_type[i] = 1;
    }
    cell_no = 1;
    ele_no = 1280;
    cell_type[0] = 1;
    ele_per_cell[0] = 1280;

    initial_chem_paras(x, y, z, id, cell_type, ele_type, cell_no, ele_per_cell);

	// Do the OpenCL calculation
    runCL(x, y, z, id, cell_type, ele_type, dt, max_ele_no, ele_no, cell_no, ele_per_cell);

	// Free up memory
	free(x);
	free(y);
	free(z);
	free(id);
    free(cell_type);
	free(ele_type);
	free(ele_per_cell);

    return 0;
}
