
// This is the kernel to run cellcontact
__kernel void cellcontact( 
__global float * xc, 
__global float * yc, 
__global float * zc, 
__const float lx, 
__const float ly, 
__const float lz, 
__global int * flag_contact,
__const int cell_no,
__global int * id_is_out) 
{
    int gId = get_global_id(0);
    float xFc = xc[gId];
    float yFc = yc[gId];
    float zFc = zc[gId];
    int flag = 0;
    int i;
    float r = 0.0f;
    for(i = 0; i < cell_no; i++)
    {
        r = (xc[i]-xFc)*(xc[i]-xFc) + (yc[i]-yFc)*(yc[i]-yFc) + (zc[i]-zFc)*(zc[i]-zFc);
        r = sqrt(r);
        if(r < 5.0f) 
        {
        flag = flag + 1;
        if(id_is_out[i] == 2) flag = flag;
        }
    }

    flag_contact[gId] = flag - 1;
    
}
