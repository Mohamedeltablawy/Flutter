#include "udf.h"



real h1;
real alfa1;
real hd1;
real alfad1;

DEFINE_CG_MOTION(piston,dt,vel,omega,time,dtime)
{real hdd;
real alfadd;
real iea;
real m;
real s;
real kh;
real ka;
real f;
real c=0.18;
real b;
real a;
real xa;



b=c/2;
a=-0.34;
xa=0.2;
iea=1.943;
m=10.887052;
kh=606.29467;
ka=2583.826;
s=m*xa*b;
f=m*iea-s*s;

Domain *d = Get_Domain(1);
Thread *t_object = Lookup_Thread(d,6);
real force[ND_ND], moment[ND_ND], cg[ND_ND]; /*initialise*/
cg[0]=b+a*b;
cg[1]=0;
cg[2]=0;
Compute_Force_And_Moment (d,t_object,cg,force,moment,TRUE);
real force_x = force[0]; /* force components of surface "Boundary_ID" */
real force_y = force[1];
real force_z = force[2];
real moment_x = moment[0]; /* moment components */
real moment_y = moment[1];
real moment_z = moment[2];

if(time<0.0055)
{h1=0;
alfa1=0;
hd1=3;
alfad1=0;
omega[2]=0;
vel[1]=hd1;
}
else
{
hdd=-iea*kh/f*h1+s*ka/f*alfa1+(iea*force_y+s*moment_z)/f;
alfadd=s*kh/f*h1-m*ka/f*alfa1+(-s*force_y-moment_z*m)/f;
hd1=hd1+hdd*dtime;
alfad1=alfad1+alfadd*dtime;
alfa1=alfa1+alfad1*dtime;
h1=h1+hd1*dtime;
omega[2]=alfad1;
vel[1]=hd1;
}
Message("dtime=%g \n",dtime);
Message("hd=%g \n",hd1);
Message("alfa=%g \n",alfa1);
Message("h=%g \n",h1);
Message("force_y=%g \n",force_y);
Message("moment_z=%g \n",moment_z);
}



DEFINE_CG_MOTION(airfole,dt,vel,omega,time,dtime)
{
if(time<0.0055)
{
omega[2]=0;
vel[1]=hd1;
}
else
{
omega[2]=alfad1;
vel[1]=hd1;
}

Message("h=%g \n",h1);
}
