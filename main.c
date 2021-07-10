
#include <calc_distance.h>
#include <stdlib.h>
#include <sys/time.h>
model m;

double get_time(void) {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}

v3 pG[10000];
int main(int argc, char *argv[])
{
    char dir[1024],file_mod[1024];
    sprintf(dir,"%s",__FILE__);
    dir[strlen(dir)-6]=0; //hack to get source directory

    //sprintf(file_mod,"%sdodeca_half_a.stl",dir);

    sprintf(file_mod,"%steapot.stl",dir);
   // sprintf(file_mod,"%s/stl_files/teapot.stl",dir);
    printf("loading %s \n", file_mod);

   modelLoad(file_mod,&m);//file_mod);

    v3 point;
    point.r[0]=-1.9;
    point.r[1]=0.0;
    point.r[2]=28.0;
    double dist=modelDistP(&m,&point);
    printf("test distance1= %le 3.495073e+000 \n", dist);

    point.r[0]=-1.9;
    point.r[1]=0.0;
    point.r[2]=28.0;
     dist=modelDistP_fast(&m,&point);
    printf("fast distance1= %le 3.495073e+000 \n", dist);



    point.r[0]=-1.9;
    point.r[1]=40.0;
    point.r[2]=28.0;
   dist=modelDistP(&m,&point);
  printf("test distance1= %le -7.663049e+000 \n", dist);

  point.r[0]=-1.9;
  point.r[1]=40.0;
  point.r[2]=28.0;
 dist=modelDistP_fast(&m,&point);
printf("fast distance1= %le -7.663049e+000 \n", dist);


  dist=0.0;
  for (int i=0;i<10000;i++)
  {
      point.r[0]=m.m_xMin +(rand()*(m.m_xMax-m.m_xMin)/RAND_MAX);
      point.r[1]=m.m_yMin +(rand()*(m.m_yMax-m.m_yMin)/RAND_MAX);
      point.r[2]=m.m_zMin +(rand()*(m.m_zMax-m.m_zMin)/RAND_MAX);

      pG[i]=point;
      }
  double t0=get_time();
  for (int i=0;i<10000;i++)
  {
    dist+=modelDistP(&m,&(pG[i]));
  }

double t1=get_time();
  printf("dit=%f  time= %e s\n",dist,t1-t0);

  dist=0.0;
   t0=get_time();
  for (int i=0;i<10000;i++)
  {
    dist+=modelDistP_fast(&m,&(pG[i]));
  }

 t1=get_time();
  printf("dit=%f  fast time= %e s\n",dist,t1-t0);

  ////for (int i=0;i<10;i++)
  //printf("freqs[%d] = %d \n", i,freqs[i]);

    return 0;
}
