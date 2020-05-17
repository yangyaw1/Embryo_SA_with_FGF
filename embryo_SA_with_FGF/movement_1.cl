// 3d distance function
float distance3(float x1, float y1, float z1, float x2, float y2, float z2, float lx, float ly, float lz)
{
    float dx = fabs(x1 - x2);
    float dy = fabs(y1 - y2);
    float dz = fabs(z1 - z2);

    return sqrt(dx*dx + dy*dy + dz*dz) + 0.0001;
}


//the intracellular potential, the morse potential in this case
float intra_potential(float r, float x, float y, float z, float * x_F1, float * y_F1, float * z_F1, int i_cycle) 
{
    if(r <= 0.0f)  return 0.0f;
    if(i_cycle <= 2 && r >= 5.0f) return 0.0f;

    float rb = 2.5f/(1.0f + 0.1*i_cycle); // balance point
    if(i_cycle <= 3) rb = 2.5f/(1.0f + 0.1*(i_cycle+1));
    if(i_cycle == 7) rb = 2.5f/(1.0f + 0.1*6);

    float r2 = 1.0f/(r*r);
    float rm = rb/r; // 1.8
    float rm6 = rm*rm*rm*rm*rm*rm;
    float m = 1.5f*(rm6*rm6*r2 - rm6*r2);
    m -= .01*r;

    if(m > 10.0) return 10.0f;
    else return 1.0f*m; 
}

//the intercellular potential, the repulsive force morse potential in this case
float inter_potential(float r, float x, float y, float z, float * x_F1, float * y_F1, float * z_F1, int i_cycle)
{
    if(r <= 0.0f || r >= 15.0f) return 0.0;

    float rb = 4.0f/(1.0f + 0.05*i_cycle);

    float r2 = 1.0f/(r*r);
    float rm = rb/r; // 4.5
    float rm6 = rm*rm*rm*rm*rm*rm;
    float m = 0.3f*(rm6*rm6*r2 - rm6*r2); // 0.1

    if(m > 10.0) return 10.0f;
    else return 1.0f*m; 
}

//the external potential
void external_potential(float x, float y, float z, float dt,
                        float * x_F1, float * y_F1, float * z_F1,
                        float lx, float ly, float lz, int id_is_out,
                        int i_cycle, int iter, float xc, float yc, float zc, int colevel, float xcc, float ycc, float zcc)
{
    float xx, yy, zz, rr, mm, xxx, yyy, zzz, rrr, mmm, x_change, y_change, z_change;
    float xx1, yy1, zz1, zz2, zz3, zz4, xx2, yy2, mm1, mm2, mm3, mm4;
    float ratiof;
    float rc, ro;
    xx1=0.7f*lx;
    yy1=0.7f*lx;
    zz1=0.7f*lx;
    xx2=0.5f*lx;
    yy2=0.5f*lx;
    mm1=0.0f;
    mm2=0.0f;
    mm3=0.0f;
    mm4=0.0f;

    // if(rr <= 1.75f) return;
    rc = distance3(xc, yc, zc, 0.0f, 0.0f, 0.0f, lx, ly, lz);
    ro = distance3(x, y, z, 0.0f, 0.0f, 0.0f, lx, ly, lz);
    if(id_is_out == 2 && i_cycle>=4 && colevel >= 1)
    {
        // for outer cells
        xx = 0.0f;
        yy = 0.0f;
        zz = 0.0f;//10.0f-20.0f/3.0f;
        rr = distance3(x, y, z, xx, yy, zz, lx, ly, lz);
        xxx = x*lx/rr;
        yyy = y*lx/rr;
        zzz = z*lx/rr;
        if(i_cycle > 3) mm = dt/0.1f*0.02f*i_cycle/rr;
        if(i_cycle == 6) mm=2.0*mm;
        xx1=0.5f*xx+0.86f*zz;
        xx2=0.5f*xx-0.86f*zz;
        yy1=0.5f*yy+0.86f*zz;
        yy2=0.5f*yy-0.86f*zz;
        zz1=0.5f*zz-0.86f*xx;
        zz2=0.5f*zz+0.86f*xx;
        zz3=0.5f*zz-0.86f*yy;
        zz4=0.5f*zz+0.86f*yy;
        mm1=0.0f;
        mm2=0.0f;
        mm3=0.0f;
        mm4=0.0f;

        mmm = 0.0f;
        if(i_cycle == 6 && rr < lx && iter<65000) mmm = (iter-40000.0f)/25000.0f*dt/0.1f*0.02f*200.0f/(lx-rr);
        if(rr < lx && iter>=65000) mmm = dt/0.1f*0.02f*200/(lx-rr);
    }
    else
    {
        // for inner cells, receive a replusive force from cavity
        xx = xcc*lx;
        yy = ycc*lx;
        zz = zcc*lx; //10.0f-20.0f/3.0f;
        rr = distance3(x, y, z, xx, yy, zz, lx, ly, lz);

        mm = dt/0.1f*0.03f*(i_cycle-3)/rr;
        xx1=0.5f*xx+0.86f*zz;
        xx2=0.5f*xx-0.86f*zz;
        yy1=0.5f*yy+0.86f*zz;
        yy2=0.5f*yy-0.86f*zz;
        zz1=0.5f*zz-0.86f*xx;
        zz2=0.5f*zz+0.86f*xx;
        zz3=0.5f*zz-0.86f*yy;
        zz4=0.5f*zz+0.86f*yy;

        rr = distance3(x, y, z, xx1, yy, zz1, lx, ly, lz);

        mm1 = dt/0.1f*0.03f*(i_cycle-3)/rr;

        rr = distance3(x, y, z, xx2, yy, zz2, lx, ly, lz);

        mm2 = dt/0.1f*0.03f*(i_cycle-3)/rr;

        rr = distance3(x, y, z, xx, yy1, zz3, lx, ly, lz);

        mm3 = dt/0.1f*0.03f*(i_cycle-3)/rr;

        rr = distance3(x, y, z, xx, yy2, zz4, lx, ly, lz);

        mm4 = dt/0.1f*0.03f*(i_cycle-3)/rr;

  
      
        // detach from boundary
        xxx = 0.0f;
        yyy = 0.0f;
        zzz = 0.0f;
        rrr = distance3(x, y, z, xxx, yyy, zzz, lx, ly, lz);
        if(lx-rrr < 0.1f && lx-rrr>0.0f && iter>50000) 
        {
        mmm = 0*dt/0.1f*0.15f*(1/(lx-rrr+0.01f)-1/(0.1f+0.01));
        if(rrr > lx && iter>50000) mmm = -0*50*2*dt/0.1f*0.15f*(1/(lx-rrr+0.01f)-1/(0.1f+0.01));
        }
        else
        {
        mmm = 0.0f; 
        }
    }
    if(mm > 50.0f) mm = 50.0f;
    if(mm1 > 50.0f) mm1 = 50.0f;
    if(mm2 > 50.0f) mm2 = 50.0f;
    if(mm3 > 50.0f) mm3 = 50.0f;
    if(mm4 > 50.0f) mm4 = 50.0f;
    if(mmm > 50.0f) mmm = 50.0f;
    if(mmm < -50.0f) mmm = -50.0f;
    if(i_cycle<6) mmm = 0.0f;
    if(i_cycle<6)
    {
    mm1=0.0f;
    mm2=0.0f;
    mm3=0.0f;
    mm4=0.0f;
    }
    if(i_cycle>=6)
    {
    ratiof=1/40;
    if(id_is_out != 2) mm=mm*ratiof;
    if(iter<=65000)
    {
    mm1=0.0f;
    mm2=0.0f;
    mm3=0.0f;
    mm4=0.0f;
    }
    if(iter>65000)
    {
    mm1=mm1*ratiof;
    mm2=mm2*ratiof;
    mm3=mm3*ratiof;
    mm4=mm4*ratiof;
    }
    }
    x_change=(x-xx)*(mm+mm3+mm4)+(x-xx1)*mm1+(x-xx2)*mm2-(x-xxx)*mmm;;
    y_change=(y-yy)*(mm+mm1+mm2)+(y-yy1)*mm3+(y-yy2)*mm4-(y-yyy)*mmm;
    z_change=(z-zz)*mm+(z-zz1)*mm1+(z-zz2)*mm2+(z-zz3)*mm3+(z-zz4)*mm4-(z-zzz)*mmm;
    *x_F1 += x_change;
    *y_F1 += y_change;
    *z_F1 += z_change;
    return;

}


// the first kernel for the RK2 method
__kernel void movement(
__global int * id, 
__global float * x, 
__global float * y, 
__global float * z, 
__global int * type, 
__const int cell_no,
__const int ele_no,
__global int * ele_per_cell, 
__global float * xc,
__global float * yc,
__global float * zc,
__const float dt, 
__global float * x_F, 
__global float * y_F, 
__global float * z_F,
__const float lx, 
__const float ly, 
__const float lz,
__const int max_ele,
__global int * cell_type,
__const int i_cycle,
__global int * id_is_out,
__const int iter,
__const int flag_division,
__const int id1,
__const int id2,
__global float * randx,
__global float * randy,
__global float * randz,
__global int * colevel,
__const float xcc,
__const float ycc,
__const float zcc,
__const int sa_on_time)
{
    int gId = get_global_id(0); //get the global ID of this work unit.
    
    int cellId = id[gId]-1;
    int icellId;

    int i, j, k, cid;
    float x_F1 = 0.0f;
    float y_F1 = 0.0f;
    float z_F1 = 0.0f;
    float cutdis;

    //create a variable to store distance
    float r = 0.0f;
    float V = 0.0f; //a temporary variable to store potentials

    // external force for empty room
    if(flag_division == 0 && i_cycle >= 4)
    {
        external_potential(x[gId], y[gId], z[gId], dt, &x_F1, &y_F1, &z_F1, lx, ly, lz, id_is_out[cellId], i_cycle, iter, xc[cellId], yc[cellId], zc[cellId], colevel[cellId], xcc, ycc, zcc);
    }


    for(k = 0; k < ele_no; k++)
    {
        i = k;

        r = distance3(x[gId], y[gId], z[gId], x[i], y[i], z[i], lx, ly, lz);
        {
           if(id[gId] == id[i])
           {
              {
                  V = intra_potential(r, x[gId], y[gId], z[gId], &x_F1, &y_F1, &z_F1, i_cycle);
              }
           }
           else
           {
              {
                  V = inter_potential(r, x[gId], y[gId], z[gId], &x_F1, &y_F1, &z_F1, i_cycle);
                  cutdis=0.0f;
                  icellId = id[i]-1;
                  if(ele_per_cell[cellId]+ele_per_cell[icellId]==20) cutdis=1.7*3.79f;
                  if(ele_per_cell[cellId]+ele_per_cell[icellId]==30) cutdis=1.7*4.189f;
                  if(ele_per_cell[cellId]+ele_per_cell[icellId]==40) cutdis=1.7*4.588f;
                  if(ele_per_cell[cellId]+ele_per_cell[icellId]==50) cutdis=1.7*4.9879f;
                  if(ele_per_cell[cellId]+ele_per_cell[icellId]==60) cutdis=1.7*5.387f;
                  if(ele_per_cell[cellId]+ele_per_cell[icellId]==80) cutdis=1.7*6.1859f;
                  if(r>=cutdis && ele_per_cell[cellId]<=40) V=0;
                  if(r<cutdis && iter >= 65000)
                  { 
                  icellId = id[i]-1;
                  int id1 = id_is_out[cellId];
                  int id2 = id_is_out[icellId];
                  if(iter >= 260000) // adhesion never occurs
                  {
                  if(id1 == 1 || id1 == 3) id1 = 1;
                  if(id2 == 1 || id2 == 3) id2 = 1;
                  }
                  if(iter < 260000) // adhesion never occurs
                  {
                  if(id1 == 0 || id1 == 1 || id1 == 3 || id1 == 4) id1 = 0;
                  if(id2 == 0 || id2 == 1 || id2 == 3 || id2 == 4) id2 = 0;
                  }
                      if(id1+id2 == 2 && id1 == id2)
                      {
                          if(V>0) V=0.25*V;
                          if(V<0) V=4*V;
                      }
                      if(id1+id2 == 0)
                      {
                          if(V>0) V=0.25*V;
                          if(V<0) V=4*V;
                      }
                      if(id1+id2 == 1)
                      {
                          if(V>0) V=0.25*V;
                          if(V<0) V=4*V;
                      }
                      if(id1+id2 == 8)
                      {
                          if(V>0) V=1.0f*V;
                          if(V<0) V=0.75f*V;
                      }
                      
                      if( id1+id2 == 4 && id1 != id2)
                      {
                          if(V>0) V=2.0f*V;
                          if(V<0) V=1.0f*V;
                      }
                      if( id1+id2 == 5 && id1 != id2)
                      {
                          if(V>0) V=1.0f*V;
                          if(V<0) V=1.0f*V;
                      }
                      if( id1+id2 == 2 && id1 != id2)
                      {
                          if(V>0) V=1*V;
                          if(V<0) V=1*V;
                      }
                      if(id1+id2 == 6)
                      {
                          if(V>0) V=6*V;
                          if(V<0) V=0.1*V;
                      }
                      if(id1+id2 == 4 && id1 == id2)
                      {
                          if(V>0) V=1*V;
                          if(V<0) V=1*V;
                      }
                  }
                  if(iter >= 50000 && iter<65000)
                  { 
                  icellId = id[i]-1;
                  if(id_is_out[cellId]+id_is_out[icellId] == 4 && id_is_out[cellId] == id_is_out[icellId])  
                  {
                  if(V>0) V=1*V;
                  if(V<0) V=1*V;
                  } 
                  }
                 
               }
           }
        }
        x_F1 += V * (x[gId] - x[i]);
        y_F1 += V * (y[gId] - y[i]);
        z_F1 += V * (z[gId] - z[i]);
    }
    x_F[gId] = x[gId]; 
    y_F[gId] = y[gId];
    z_F[gId] = z[gId];
    if(i_cycle != 7 || id_is_out[cellId] != 2)
    {
    if(i_cycle>=4)
    {
    x_F[gId] = .5 * dt * (x_F1 + 0.25*randx[gId]*x_F1) + x[gId]; //save the calculated values
    y_F[gId] = .5 * dt * (y_F1 + 0.25*randy[gId]*y_F1) + y[gId];
    z_F[gId] = .5 * dt * (z_F1 + 0.25*randz[gId]*z_F1) + z[gId];
    }
    if(i_cycle<4)
    {
    x_F[gId] = .5 * dt * (x_F1 + 0.25*randx[gId]*x_F1) + x[gId]; //save the calculated values
    y_F[gId] = .5 * dt * (y_F1 + 0.25*randy[gId]*y_F1) + y[gId];
    z_F[gId] = .5 * dt * (z_F1 + 0.25*randz[gId]*z_F1) + z[gId];
    }
    if(flag_division == 1)
    {
        x_F[gId] = .1 * dt * x_F1 + x[gId]; //save the calculated values
        y_F[gId] = .1 * dt * y_F1 + y[gId];
        z_F[gId] = .1 * dt * z_F1 + z[gId];
    }

    r = distance3(x_F[gId], y_F[gId], z_F[gId], 0.0f, 0.0f, 0.0f, lx, ly, lz);
    float RR = lx;
    if(r > RR)
    {
        x_F[gId] = x_F[gId]*RR/r;
        y_F[gId] = y_F[gId]*RR/r;
        z_F[gId] = z_F[gId]*RR/r;
    }
    else if(r > RR-1.5f && id_is_out[cellId]!=2 && iter>2275000)
    {
        x_F[gId] = x_F[gId]*(RR-1.5)/r;
        y_F[gId] = y_F[gId]*(RR-1.5)/r;
        z_F[gId] = z_F[gId]*(RR-1.5)/r;
    }
    }

}

