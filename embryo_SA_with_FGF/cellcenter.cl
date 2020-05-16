
// This is the kernel to run cellcenter
__kernel void cellcenter( 
__global int * id, 
__global float * x, 
__global float * y, 
__global float * z, 
__global int * ele_per_cell, 
__global float * xc,
__global float * yc,
__global float * zc,
__const int max_ele,
__const float lx,
__const float ly,
__const float lz,
__const int ele_no,
__global int * flag_outer,
__const int iter) 
{
    int gId = get_global_id(0); //get the global ID of this work unit. This will be the cell number of interest
    
    float cx = 0.0f; //Calculate the center of the top of the stem cell in the x and y direction 
    float cy = 0.0f;
    float cz = 0.0f;
    int count_tmp = 0;
    int i, k, kk;
    int flag = 0;
    float r;

    for(i = 0; i < ele_no; i++)
    {
        if((id[i]-1) == gId)
        {
            cx += x[i];
            cy += y[i];
            cz += z[i];
            count_tmp++;
            // if(flag == 0)
            {
                r = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
                r = sqrt(r);
                if(lx-r < 1.25f) 
                {
                    flag = flag + 1;
                }
            }
        }
    }

    flag_outer[gId] = flag;
    if(count_tmp == ele_per_cell[gId])
    {
        cx = cx/count_tmp;
        cy = cy/count_tmp;
        cz = cz/count_tmp;
        xc[gId] = cx;
        yc[gId] = cy;
        zc[gId] = cz;
    }

}
