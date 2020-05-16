      rr = distance3(x, y, z, xx1, yy, zz1, lx, ly, lz);

        mm1 = dt/0.1f*0.03f*(i_cycle-3)/rr;

        rr = distance3(x, y, z, -xx1, yy, zz1, lx, ly, lz);

        mm2 = dt/0.1f*0.03f*(i_cycle-3)/rr;

        rr = distance3(x, y, z, xx, yy1, zz1, lx, ly, lz);

        mm3 = dt/0.1f*0.03f*(i_cycle-3)/rr;

        rr = distance3(x, y, z, xx, -yy1, zz1, lx, ly, lz);

        mm4 = dt/0.1f*0.03f*(i_cycle-3)/rr; 


    x_change=(x-xx)*(mm+mm3+mm4)+(x—xx1)*mm1+(x+xx1)*mm2-(x-xxx)*mmm;
    y_change=(y-yy)*(mm+mm1+mm2)+(y-yy1)*mm3+(y+yy1)*mm4-(y-yyy)*mmm;
    z_change=(z-zz)*mm+(z-zz1)*(mm1+mm2+mm3+mm4)-(z-zzz)*mmm;