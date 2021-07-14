#ifndef CALC_DISTANCE_H
#define CALC_DISTANCE_H

#include <stdio.h>
#include <math.h>
#include <x86intrin.h>
#define N_TRIS_MAX 20000

typedef struct{
    double r[3];//, r[1], r[2];
} v3;

typedef struct{
    v3 m_p[3];	//coordinates of three points

    v3 shP[3];	//shifted points - coordinates after transformation
    double m_vert_angle[3]; //angles at the vertices (in rad)
    v3 normal; //normal vector
    v3 center; //point of triangle centroid

    double matrix[4][4]; //transformation matrix (translation + rotation)
} tri;

typedef struct {
    tri m_tris[N_TRIS_MAX];
    int n_m_tris;

    double m_xMin, m_xMax, m_yMin, m_yMax, m_zMin, m_zMax;

    double w_x0, w_y0, w_z0, w_x1, w_y1, w_z1;
    double *matr_1; //line1
    double *p_tr; //transformed points
    double *shP_1;
    double *shP_2;
    double *ang;
} model;



v3 v3Init(char* facet) {
    v3 vec;
    char f1[4] = { facet[0],
        facet[1], facet[2], facet[3] };
    char f2[4] = { facet[4],
        facet[5], facet[6], facet[7] };
    char f3[4] = { facet[8],
        facet[9], facet[10], facet[11] };

    float xx = *((float*)f1);
    float yy = *((float*)f2);
    float zz = *((float*)f3);

    vec.r[0] = (double)(xx);
    vec.r[1] = (double)(yy);
    vec.r[2] = (double)(zz);
    return vec;
}

v3 multiplyMatrPoint(double A[4][4], v3 * b) {
    v3 res;

    res.r[0] = A[0][0]*b->r[0] + A[0][1]*b->r[1] + A[0][2]*b->r[2] + A[0][3];
    res.r[1] = A[1][0]*b->r[0] + A[1][1]*b->r[1] + A[1][2]*b->r[2] + A[1][3];
    res.r[2] = A[2][0]*b->r[0] + A[2][1]*b->r[1] + A[2][2]*b->r[2] + A[2][3];

    return res;
}


void multiplyMatrixes(double A[4][4], double B[4][4]) {
    double temp[4][4];

    for(int i=0;i<4;i++) {
        for(int j=0;j<4;j++) {
            temp[i][j] = 0;
            for(int m=0;m<4;m++) {
                temp[i][j] += A[i][m]*B[m][j];
            }
        }
    }

    for(int i=0;i<4;i++) {
        for(int j=0;j<4;j++) {
            B[i][j] = temp[i][j];
        }
    }
}

double len(v3 v) {
    return sqrt(v.r[0]*v.r[0] + v.r[1]*v.r[1] + v.r[2]*v.r[2]);
}

double dotProd(v3 v1, v3 v2) {
    return v1.r[0]*v2.r[0] + v1.r[1]*v2.r[1] + v1.r[2]*v2.r[2];
}

#define DOT(v1,v2) (v1.r[0]*v2.r[0] + v1.r[1]*v2.r[1] + v1.r[2]*v2.r[2])
#define LEN2(v) (v.r[0]*v.r[0] + v.r[1]*v.r[1] + v.r[2]*v.r[2])

tri triInit(v3 p1, v3 p2, v3 p3) {
    tri t;
    t.m_p[0] = p1;
    t.m_p[1] = p2;
    t.m_p[2] = p3;

    //cross product [p1p2 x p1p3]
    t.normal.r[0] = (p2.r[1] - p1.r[1])*(p3.r[2] - p1.r[2]) - (p2.r[2] - p1.r[2])*(p3.r[1] - p1.r[1]);
    t.normal.r[1] = (p2.r[2] - p1.r[2])*(p3.r[0] - p1.r[0]) - (p2.r[0] - p1.r[0])*(p3.r[2] - p1.r[2]);
    t.normal.r[2] = (p2.r[0] - p1.r[0])*(p3.r[1] - p1.r[1]) - (p2.r[1] - p1.r[1])*(p3.r[0] - p1.r[0]);

    double l = len(t.normal);
    t.normal.r[0] /= l;
    t.normal.r[1] /= l;
    t.normal.r[2] /= l;

    t.center.r[0] = (p1.r[0] + p2.r[0] + p3.r[0]) / 3.;
    t.center.r[1] = (p1.r[1] + p2.r[1] + p3.r[1]) / 3.;
    t.center.r[2] = (p1.r[2] + p2.r[2] + p3.r[2]) / 3.;

    for(int i=0;i<4;i++) {
        for(int j=0;j<4;j++) {
            if (i != j) t.matrix[i][j] = 0.;
            else t.matrix[i][j] = 1.;
        }
    }

    //translation for p1 is in the origin
    t.matrix[0][3] = -p1.r[0];
    t.matrix[1][3] = -p1.r[1];
    t.matrix[2][3] = -p1.r[2];

    for (int i = 0; i < 3; i++) {
        t.shP[i].r[0] = t.m_p[i].r[0] - p1.r[0];
        t.shP[i].r[1] = t.m_p[i].r[1] - p1.r[1];
        t.shP[i].r[2] = t.m_p[i].r[2] - p1.r[2];
    }

    double rot_matrix[4][4];   //rotation matrix

    for(int i=0;i<4;i++) {
        for(int j=0;j<4;j++) {
            if(i!=j) rot_matrix[i][j] = 0.;
            else rot_matrix[i][j] = 1.;
        }
    }

    //angle between vector p1p2 and plane y=0
    double alpha = acos(t.shP[1].r[1] / sqrt(t.shP[1].r[0]*t.shP[1].r[0] + t.shP[1].r[1]*t.shP[1].r[1]));
    if (t.shP[1].r[0] > 0) alpha = -alpha;
    //cout << "alpha=" << alpha*180/M_PI << endl;
    //rotation for p1p2 is on the xz-plane
    rot_matrix[0][0] = cos(-alpha + M_PI/2.);
    rot_matrix[0][1] = -sin(-alpha + M_PI/2.);
    rot_matrix[1][0] = sin(-alpha + M_PI/2.);
    rot_matrix[1][1] = cos(-alpha + M_PI/2.);

    //multiplyMatrVect(rot_matrix, t.shP);  //should be p2.y=0
    for (int i = 0; i < 3; i++) t.shP[i] = multiplyMatrPoint(rot_matrix, &(t.shP[i]));

    multiplyMatrixes(rot_matrix, t.matrix);

    for(int i=0;i<4;i++) {
        for(int j=0;j<4;j++) {
            if(i!=j) rot_matrix[i][j] = 0.;
            else rot_matrix[i][j] = 1.;
        }
    }

    //angle between vector p1p2 and plane x=0
    double beta = acos(t.shP[1].r[0] / sqrt(t.shP[1].r[0]*t.shP[1].r[0] + t.shP[1].r[2]*t.shP[1].r[2]));
    if (t.shP[1].r[2] < 0) beta = -beta;

    //rotation for p1p2 and z-axis are co-directed
    rot_matrix[0][0] = cos(beta - M_PI/2.);
    rot_matrix[0][2] = sin(beta - M_PI/2.);
    rot_matrix[2][0] = -sin(beta - M_PI/2.);
    rot_matrix[2][2] = cos(beta - M_PI/2.);

    //multiplyMatrVect(rot_matrix, t.shP);  //should be p2.x=0
    for (int i = 0; i < 3; i++) t.shP[i] = multiplyMatrPoint(rot_matrix, &(t.shP[i]));

    multiplyMatrixes(rot_matrix, t.matrix);

    for(int i=0;i<4;i++) {
        for(int j=0;j<4;j++) {
            if(i!=j) rot_matrix[i][j] = 0.;
            else rot_matrix[i][j] = 1.;
        }
    }

    //angle between vector p1p3 and plane x=0
    double gamma = acos(t.shP[2].r[0] / sqrt(t.shP[2].r[0]*t.shP[2].r[0] + t.shP[2].r[1]*t.shP[2].r[1]));
    if (t.shP[2].r[1] < 0) gamma = -gamma;

    //rotation for p3p1 is on the yz-plane
    rot_matrix[0][0] = cos(-gamma + M_PI/2.);
    rot_matrix[0][1] = -sin(-gamma + M_PI/2.);
    rot_matrix[1][0] = sin(-gamma + M_PI/2.);
    rot_matrix[1][1] = cos(-gamma + M_PI/2.);

    //multiplyMatrVect(rot_matrix, t.shP); //should be p3.x=0
    for (int i = 0; i < 3; i++) t.shP[i] = multiplyMatrPoint(rot_matrix, &(t.shP[i]));

    multiplyMatrixes(rot_matrix, t.matrix);

  //find incident angles at vertices;
    v3 a,b;
     a.r[0]=t.m_p[1].r[0] - t.m_p[0].r[0];
     a.r[1]=t.m_p[1].r[1] - t.m_p[0].r[1];
     a.r[2]=t.m_p[1].r[2] - t.m_p[0].r[2];

     b.r[0]=t.m_p[2].r[0] - t.m_p[0].r[0];
     b.r[1]=t.m_p[2].r[1] - t.m_p[0].r[1];
     b.r[2]=t.m_p[2].r[2] - t.m_p[0].r[2];
    t.m_vert_angle[0]=acos( dotProd(a,b) / (len(a)*len(b)) );

     a.r[0]=t.m_p[0].r[0] - t.m_p[1].r[0];
     a.r[1]=t.m_p[0].r[1] - t.m_p[1].r[1];
     a.r[2]=t.m_p[0].r[2] - t.m_p[1].r[2];

     b.r[0]=t.m_p[2].r[0] - t.m_p[1].r[0];
     b.r[1]=t.m_p[2].r[1] - t.m_p[1].r[1];
     b.r[2]=t.m_p[2].r[2] - t.m_p[1].r[2];

    t.m_vert_angle[1]=acos( dotProd(a,b) / (len(a)*len(b)) );

     a.r[0]=t.m_p[1].r[0] - t.m_p[2].r[0];
     a.r[1]=t.m_p[1].r[1] - t.m_p[2].r[1];
     a.r[2]=t.m_p[1].r[2] - t.m_p[2].r[2];

     b.r[0]=t.m_p[0].r[0] - t.m_p[2].r[0];
     b.r[1]=t.m_p[0].r[1] - t.m_p[2].r[1];
     b.r[2]=t.m_p[0].r[2] - t.m_p[2].r[2];

    t.m_vert_angle[2]=acos( dotProd(a,b) / (len(a)*len(b)) );

    if (fabs(t.m_vert_angle[0] + t.m_vert_angle[1] + t.m_vert_angle[2]- 3.1415926535)>1e-5 )
        printf("WARNING triangle angles sum is %f \n",t.m_vert_angle[0] + t.m_vert_angle[1] +t.m_vert_angle[2]);


    return t;
}

double distP_tri(tri* t, v3* point, double* dot_out) {
    double dist = -1;
    v3 p = multiplyMatrPoint(t->matrix, point);

    double E12, E23, E31; //edges equation results: E<0 -> point on left side, E>0 -> point on right side
    //E(y,z) = (y-Y)*dZ - (z-Z)*dY
    (*dot_out)=p.r[0];
    E12 = p.r[1]*t->shP[1].r[2];
    E23 = p.r[1]*(t->shP[2].r[2] - t->shP[1].r[2]) - (p.r[2] - t->shP[1].r[2])*t->shP[2].r[1];
    E31 = p.r[1]*(-t->shP[2].r[2]) - p.r[2]*(-t->shP[2].r[1]);
//printf("E12=%f E23=%f E31=%f ------ \n", E12,E23,E31);

    if(E12>-1.e-15 && E23>-1.e-15 && E31>-1e-15) {   //inside
        dist = fabs(p.r[0]);
        return dist;
    }
    double mult = 0;
    if(E12 < 0) { //left side of p1p2
        mult = dotProd(t->shP[1], p);
        if(mult < 0) {  //p1 is closest
            (*dot_out)*=t->m_vert_angle[0];
            return len(p);
        }
        else {
            mult /= fabs(t->shP[1].r[2]);
            if (mult <= fabs(t->shP[1].r[2])) {     //p1p2 is closest
                return sqrt(p.r[0]*p.r[0] + p.r[1]*p.r[1]);
            }
            else {  //p2 is closest
                v3 temp;
                (*dot_out)*=t->m_vert_angle[1];
                temp.r[0] = p.r[0] - t->shP[1].r[0];
                temp.r[1] = p.r[1] - t->shP[1].r[1];
                temp.r[2] = p.r[2] - t->shP[1].r[2];

                return len(temp);
            }
        }
    }

    if(E31 < 0) {  //right side p1p2, left side p3p1
        mult = dotProd(t->shP[2], p);
        if(mult < 0) {  //p1 is closest
            (*dot_out)*=t->m_vert_angle[0];
            return len(p);
        }
        else {
            mult /= fabs(len(t->shP[2]));
            if (mult <= fabs(len(t->shP[2]))) {     //p1p3 is closest
                double c2 = p.r[1]*p.r[1] + p.r[2]*p.r[2] - mult*mult;
                return sqrt(p.r[0]*p.r[0] + c2);
            }
            else {  //p3 is closest
                v3 temp;
                (*dot_out)*=t->m_vert_angle[2];
                temp.r[0] = p.r[0] - t->shP[2].r[0];
                temp.r[1] = p.r[1] - t->shP[2].r[1];
                temp.r[2] = p.r[2] - t->shP[2].r[2];

                return len(temp);
            }
        }
    }

    if(E23 < 0) {  //left side p2p3, right side of other
        v3 temp1, temp2;

        temp1.r[0] = t->shP[2].r[0] - t->shP[1].r[0];
        temp1.r[1] = t->shP[2].r[1] - t->shP[1].r[1];
        temp1.r[2] = t->shP[2].r[2] - t->shP[1].r[2];

        temp2.r[0] = p.r[0] - t->shP[1].r[0];
        temp2.r[1] = p.r[1] - t->shP[1].r[1];
        temp2.r[2] = p.r[2] - t->shP[1].r[2];

        mult = dotProd(temp1, temp2);
        if(mult < 0) {  //p2 is closest
            (*dot_out)*=t->m_vert_angle[1];
            return len(temp2);
        }
        else {
            mult /= len(temp1);
            if (mult <= len(temp1)) {     //p2p3 is closest
                double c2 = pow(temp2.r[1], 2.) + pow(temp2.r[2], 2) - mult*mult;
                return sqrt(p.r[0]*p.r[0] + c2);
            }
            else {  //p3 is closest
                v3 temp;
                (*dot_out)*=t->m_vert_angle[2];
                temp.r[0] = p.r[0] - t->shP[2].r[0];
                temp.r[1] = p.r[1] - t->shP[2].r[1];
                temp.r[2] = p.r[2] - t->shP[2].r[2];

                return len(temp);
            }
        }
    }
    return dist;
}


double distP_tri_optimized(tri* t, v3* point, double* dot_out) {
    double dist = -1;

    v3 p;
  /*  p.r[0] = t->matrix[0][0]*point->r[0] + t->matrix[0][1]*point->r[1] + t->matrix[0][2]*point->r[2] + t->matrix[0][3];
    p.r[1] = t->matrix[1][0]*point->r[0] + t->matrix[1][1]*point->r[1] + t->matrix[1][2]*point->r[2] + t->matrix[1][3];
    p.r[2] = t->matrix[2][0]*point->r[0] + t->matrix[2][1]*point->r[1] + t->matrix[2][2]*point->r[2] + t->matrix[2][3];
*/
    double * ai=&(t->matrix[0][0]);
    double * ai1=&(t->matrix[1][0]);
    double * ai2=&(t->matrix[2][0]);
    double *r1=&(point->r[0]);
       p.r[0] = (*ai)*(*r1) + (*(ai+1))*(*(r1+1)) + (*(ai+2))*(*(r1+2)) + (*(ai+3));
       p.r[1] = (*ai1)*(*r1) + (*(ai1+1))*(*(r1+1)) + (*(ai1+2))*(*(r1+2)) + (*(ai1+3));
       p.r[2] = (*ai2)*(*r1) + (*(ai2+1))*(*(r1+1)) + (*(ai2+2))*(*(r1+2)) + (*(ai2+3));

    double E12, E23, E31; //edges equation results: E<0 -> point on left side, E>0 -> point on right side

    (*dot_out)=p.r[0]; //normal signed distance
    //E(y,z) = (y-Y)*dZ - (z-Z)*dY

    E12 = p.r[1]*t->shP[1].r[2];
    E23 = p.r[1]*(t->shP[2].r[2] - t->shP[1].r[2]) - (p.r[2] - t->shP[1].r[2])*t->shP[2].r[1];
    E31 = p.r[1]*(-t->shP[2].r[2]) - p.r[2]*(-t->shP[2].r[1]);


    if(E12>-1.e-15 && E23>-1.e-15 && E31>-1e-15) {   //inside
        return p.r[0] * p.r[0];
    }
    double mult = 0;
    if(E12 < 0) { //left side of p1p2
        mult = DOT(t->shP[1], p);
        if(mult < 0) {  //p1 is closest
            (*dot_out)*=t->m_vert_angle[0];
            return LEN2(p);
        }
        else {
            if (mult <= (t->shP[1].r[2])*(t->shP[1].r[2])) {     //p1p2 is closest
                return (p.r[0]*p.r[0] + p.r[1]*p.r[1]);
            }
            else {  //p2 is closest
                v3 temp;
                (*dot_out)*=t->m_vert_angle[1];
                temp.r[0] = p.r[0] - t->shP[1].r[0];
                temp.r[1] = p.r[1] - t->shP[1].r[1];
                temp.r[2] = p.r[2] - t->shP[1].r[2];
                return LEN2(temp);
            }
        }
    }

    if(E31 < 0) {  //right side p1p2, left side p3p1
        mult = DOT(t->shP[2], p);
        if(mult < 0) {  //p1 is closest
            (*dot_out)*=t->m_vert_angle[0];
            return LEN2(p);
        }
        else {
            double l2=LEN2(t->shP[2]);
            if (mult <= l2) {     //p1p3 is closest
                double c2 = p.r[1]*p.r[1] + p.r[2]*p.r[2] - mult*mult/l2;
                return (p.r[0]*p.r[0] + c2);
            }
            else {  //p3 is closest
                v3 temp;
                (*dot_out)*=t->m_vert_angle[2];
                temp.r[0] = p.r[0] - t->shP[2].r[0];
                temp.r[1] = p.r[1] - t->shP[2].r[1];
                temp.r[2] = p.r[2] - t->shP[2].r[2];
                return LEN2(temp);
            }
        }
    }

    if(E23 < 0) {  //left side p2p3, right side of other
        v3 temp1, temp2;

        temp1.r[0] = t->shP[2].r[0] - t->shP[1].r[0];
        temp1.r[1] = t->shP[2].r[1] - t->shP[1].r[1];
        temp1.r[2] = t->shP[2].r[2] - t->shP[1].r[2];

        temp2.r[0] = p.r[0] - t->shP[1].r[0];
        temp2.r[1] = p.r[1] - t->shP[1].r[1];
        temp2.r[2] = p.r[2] - t->shP[1].r[2];

        mult = DOT(temp1, temp2);
        if(mult < 0) {  //p2 is closest
            (*dot_out)*=t->m_vert_angle[1];
            return LEN2(temp2);
        }
        else {
            double l2= LEN2(temp1);
            if (mult <= l2) {     //p2p3 is closest
                double c2 = temp2.r[1]*temp2.r[1] + temp2.r[2]*temp2.r[2] - mult*mult/l2;
                return (p.r[0]*p.r[0] + c2);
            }
            else {  //p3 is closest
                v3 temp;
                (*dot_out)*=t->m_vert_angle[2];
                temp.r[0] = p.r[0] - t->shP[2].r[0];
                temp.r[1] = p.r[1] - t->shP[2].r[1];
                temp.r[2] = p.r[2] - t->shP[2].r[2];
                return LEN2(temp);
            }
        }
    }
    return dist;
}


void getMinMax(model* m)
{
    m->m_xMin = 1e10;    m->m_xMax = -1e10;
    m->m_yMin = 1e10;    m->m_yMax = -1e10;
    m->m_zMin = 1e10;    m->m_zMax = -1e10;

    for (int i = 0; i < m->n_m_tris; i++)
    {
        for (int j = 0; j<3; j++)
        {
            if (m->m_tris[i].m_p[j].r[0] > m->m_xMax) m->m_xMax = m->m_tris[i].m_p[j].r[0];
            if (m->m_tris[i].m_p[j].r[0] < m->m_xMin) m->m_xMin = m->m_tris[i].m_p[j].r[0];

            if (m->m_tris[i].m_p[j].r[1] > m->m_yMax) m->m_yMax = m->m_tris[i].m_p[j].r[1];
            if (m->m_tris[i].m_p[j].r[1] < m->m_yMin) m->m_yMin = m->m_tris[i].m_p[j].r[1];

            if (m->m_tris[i].m_p[j].r[2] > m->m_zMax) m->m_zMax = m->m_tris[i].m_p[j].r[2];
            if (m->m_tris[i].m_p[j].r[2] < m->m_zMin) m->m_zMin = m->m_tris[i].m_p[j].r[2];
        }
    }
    return;
}

void modelLoad(const char *fname, model* m) {

    m->n_m_tris = 0;

    FILE *myFile;
    myFile = fopen(fname, "rb");

    char header_info[80] = "";

    char nTri[4];
    unsigned long nTriLong;


    //read 80 byte header
    if (myFile) {
        for (int i = 0; i < 80;i++) header_info[i] = fgetc(myFile);
        printf("header:");
        printf("%s \n", header_info);
    }
    else{
        printf("error \n");
    }

    //read 4-byte ulong
    if (myFile) {
        //myFile.read (nTri, 4);
        for (int i = 0; i < 4; i++) nTri[i] = fgetc(myFile);
        nTriLong = *((unsigned long*)nTri) ;
        printf("n Tri:");
        printf("%d \n", nTriLong);
    }
    else{
        printf("error \n");
    }
    if ((int)(nTriLong) > N_TRIS_MAX)
    {
        printf("NTRIS %d > NTRIMAX %d \n", nTriLong, N_TRIS_MAX);
        return m;
    }

    //now read in all the triangles
    for(int i = 0; i < (int)(nTriLong); i++) {
        char facet[50];
        if (myFile) {
            //read one 50-byte triangle
            for (int j = 0; j < 50; j++) facet[j] = fgetc(myFile);

            //populate each point of the triangle
            //using v3::v3(char* bin);
            //facet + 12 skips the triangle's unit normal
            v3 p1 = v3Init(facet + 12);
            v3 p2 = v3Init(facet + 24);
            v3 p3 = v3Init(facet + 36);

            //add a new triangle to the array
            m->m_tris[m->n_m_tris] = triInit(p1, p2, p3);
            m->n_m_tris++;
        }
    }

     getMinMax(m);

    printf("model loaded name=%s tri_count=%d \n xmin=%f xmax=%f \n ymin=%f ymax=%f \n zmin=%f zmax=%f \n", fname, m->n_m_tris,
        m->m_xMin, m->m_xMax, m->m_yMin, m->m_yMax, m->m_zMin, m->m_zMax);

    double dx,dy,dz;
    m->w_x0 = m->m_xMin;
    m->w_x1 = m->m_xMax;
    m->w_y0 = m->m_yMin;
    m->w_y1 = m->m_yMax;
    m->w_z0 = m->m_zMin;
    m->w_z1 = m->m_zMax;

    dx = m->w_x1 - m->w_x0; dy = m->w_y1 - m->w_y0; dz = m->w_z1 - m->w_z0;

    m->w_x0 -= dx; m->w_x1 += dx;
    m->w_y0 -= dy; m->w_y1 += dy;
    m->w_z0 -= dz; m->w_z1 += dz;

//for fast sse:
    m->matr_1 = (double*)_mm_malloc(sizeof(double)*(m->n_m_tris*16), 32);
    m->p_tr =   (double*)_mm_malloc(sizeof(double)*(m->n_m_tris*4), 32);

    m->shP_1 =   (double*)_mm_malloc(sizeof(double)*(m->n_m_tris*4), 32);
    m->shP_2 =   (double*)_mm_malloc(sizeof(double)*(m->n_m_tris*4), 32);

    m->ang =   (double*)_mm_malloc(sizeof(double)*(m->n_m_tris*4), 32);


    for(int i = 0; i <  m->n_m_tris; i++) {
        //4x 4y 4z 4w
        m->matr_1[i*16]=m->m_tris[i].matrix[0][0];
        m->matr_1[i*16+1]=m->m_tris[i].matrix[1][0];
        m->matr_1[i*16+2]=m->m_tris[i].matrix[2][0];
        m->matr_1[i*16+3]=0.0;

        m->matr_1[i*16+4]=m->m_tris[i].matrix[0][1];
        m->matr_1[i*16+5]=m->m_tris[i].matrix[1][1];
        m->matr_1[i*16+6]=m->m_tris[i].matrix[2][1];
        m->matr_1[i*16+7]=0.0;

        m->matr_1[i*16+8]=m->m_tris[i].matrix[0][2];
        m->matr_1[i*16+9]=m->m_tris[i].matrix[1][2];
        m->matr_1[i*16+10]=m->m_tris[i].matrix[2][2];
        m->matr_1[i*16+11]=0.0;

        m->matr_1[i*16+12]=m->m_tris[i].matrix[0][3];
        m->matr_1[i*16+13]=m->m_tris[i].matrix[1][3];
        m->matr_1[i*16+14]=m->m_tris[i].matrix[2][3];
        m->matr_1[i*16+15]=0.0;

        m->p_tr[i*4]=1.0;
        m->p_tr[i*4+1]=1.0;
        m->p_tr[i*4+2]=1.0;
        m->p_tr[i*4+3]=1.0;
        //shP[1].r[0] (wall normal) is always 0 isn't it?
        //shP[1].r[0]=0 shP[1].r[1]=0 shP[2].r[0]=0
        m->shP_1[i*4]   = m->m_tris[i].shP[1].r[0];
        m->shP_1[i*4+1] = m->m_tris[i].shP[1].r[1];
        m->shP_1[i*4+2] = m->m_tris[i].shP[1].r[2];
        m->shP_1[i*4+3] = (m->m_tris[i].shP[2].r[2] -
                      m->m_tris[i].shP[1].r[2]);

        m->shP_2[i*4]   = m->m_tris[i].shP[2].r[0];
        m->shP_2[i*4+1] = m->m_tris[i].shP[2].r[1];
        m->shP_2[i*4+2] = m->m_tris[i].shP[2].r[2];
        m->shP_2[i*4+3] =m->m_tris[i].shP[1].r[2] *
                        m->m_tris[i].shP[2].r[1];


        m->ang[i*4]   = m->m_tris[i].m_vert_angle[0];
        m->ang[i*4+1] = m->m_tris[i].m_vert_angle[1];
        m->ang[i*4+2] = m->m_tris[i].m_vert_angle[2];
        m->ang[i*4+3] = 0.0;
//printf("shp1 =%f shp2 =%f shp3 =%f \n", m->m_tris[i].shP[2].r[0],m->m_tris[i].shP[2].r[1],m->m_tris[i].shP[2].r[2]);
    }

    return;
}
double modelDistP(model* m, v3* point)
{
    double min_dist =1e10;
    double sum_dot=0.0;//(0,0,0);
    double delta=1e-5;
    for (int i=0; i<m->n_m_tris;i++)
    {

        double  dot_n;
        double l=distP_tri(&(m->m_tris[i]),point,&dot_n);

        if (fabs(min_dist)>=fabs(l)-delta)
        {
            if  (fabs(l-min_dist)>=delta)
                sum_dot=dot_n;
            else
                sum_dot+=dot_n;

            min_dist=l;
        }
    }
    return min_dist*(2.0*(sum_dot>0.0)-1.0);
}

double  modelDistP_fast(model* m, v3* point)
{
    double min_dist =1e10;
    double sum_dot=0.0;//(0,0,0);
    double delta=1e-5*1e-5;
    for (int i=0; i<m->n_m_tris;i++)
    {

        double  dot_n;
        double l=distP_tri_optimized(&(m->m_tris[i]),point,&dot_n);

        if (min_dist>=l-delta)
        {
            if  (fabs(l-min_dist)>=delta)
                sum_dot=dot_n;
            else
                sum_dot+=dot_n;

            min_dist=l;
        }
    }
    return sqrt(min_dist)*(2.0*(sum_dot>0.0)-1.0);
}


void dot4x4(double *aosoa, double *b, double *out) {
    __m256d vx = _mm256_load_pd(&aosoa[0]);
    __m256d vy = _mm256_load_pd(&aosoa[4]);
    __m256d vz = _mm256_load_pd(&aosoa[8]);
    __m256d vw = _mm256_load_pd(&aosoa[12]);
    __m256d brod1 = _mm256_set1_pd(b[0]);
    __m256d brod2 = _mm256_set1_pd(b[1]);
    __m256d brod3 = _mm256_set1_pd(b[2]);
    __m256d brod4 = _mm256_set1_pd(b[3]);
    __m256d dot4 = _mm256_add_pd(
        _mm256_add_pd(_mm256_mul_pd(brod1, vx), _mm256_mul_pd(brod2, vy)),
        _mm256_add_pd(_mm256_mul_pd(brod3, vz), _mm256_mul_pd(brod4, vw)));
    _mm256_store_pd(out, dot4);

}

/*#define X_ (i*4)
#define Y_ (i*4+1)
#define Z_ (i*4+2)
#define W_ (i*4+3)
*/
double dist_sse(model* m,int i,double* dot_out)
{

        int X_=i*4;
        int Y_=i*4+1;
        int Z_=i*4+2;
        int W_=i*4+3;

  //  double E12, E23, E31; //edges equation results: E<0 -> point on left side, E>0 -> point on right side

    (*dot_out)=m->p_tr[X_]; //normal signed distance
    //E(y,z) = (y-Y)*dZ - (z-Z)*dY

    /*
    E12 = m->p_tr[X_]*t->shP[1].r[2];
    E23 = m->p_tr[Y_]*(t->shP[2].r[2] - t->shP[1].r[2]) - (m->p_tr[Z_] - t->shP[1].r[2])*t->shP[2].r[1];
    E31 = m->p_tr[Y_]*(-t->shP[2].r[2]) - m->p_tr[Z_]*(-t->shP[2].r[1]);
*/

        double aa=m->p_tr[Z_]*m->shP_2[Y_];
   double E12 = m->p_tr[Y_]*m->shP_1[Z_];
   double E23 =  m->p_tr[Y_]*m->shP_1[W_] - aa + m->shP_2[W_];
   double E31 = -m->p_tr[Y_]*m->shP_2[Z_] + aa;

   // printf("E12=%f E23=%f E31=%f \n",E12,E23,E31);

    if(E12>-1.e-15 && E23>-1.e-15 && E31>-1e-15) {   //inside
        return m->p_tr[X_] * m->p_tr[X_];

    }
    double mult = 0;
    if(E12 < 0) { //left side of p1p2
        mult = m->p_tr[Z_]*m->shP_1[Z_];//DOT(t->shP[1], p);
        if(mult < 0) {  //p1 is closest
            (*dot_out)*=m->ang[i*4];//t->m_vert_angle[0];
            return m->p_tr[X_]*m->p_tr[X_] + m->p_tr[Y_]*m->p_tr[Y_] + m->p_tr[Z_]*m->p_tr[Z_];

        }
        else {
            if (mult <= m->shP_1[Z_]*m->shP_1[Z_]) {     //p1p2 is closest
                return (m->p_tr[X_]*m->p_tr[X_] + m->p_tr[Y_]*m->p_tr[Y_]);

            }
            else {  //p2 is closest
                v3 temp;
                (*dot_out)*=m->ang[i*4+1];
                temp.r[0] = m->p_tr[X_];
                temp.r[1] = m->p_tr[Y_];
                temp.r[2] = m->p_tr[Z_] - m->shP_1[Z_];
                return  LEN2(temp);

            }
        }
    }

    if(E31 < 0) {  //right side p1p2, left side p3p1
        mult =  m->p_tr[Y_]*m->shP_2[Y_] + m->p_tr[Z_]*m->shP_2[Z_];//DOT(t->shP[2], p);
        if(mult < 0) {  //p1 is closest
            (*dot_out)*=m->ang[i*4];
            return (m->p_tr[X_]*m->p_tr[X_] + m->p_tr[Y_]*m->p_tr[Y_] + m->p_tr[Z_]*m->p_tr[Z_]);

        }
        else {
            double l2= m->shP_2[Y_]*m->shP_2[Y_] + m->shP_2[Z_]*m->shP_2[Z_];//LEN2(t->shP[2]);
            if (mult <= l2) {     //p1p3 is closest
                double c2 = m->p_tr[Y_]*m->p_tr[Y_] + m->p_tr[Z_]*m->p_tr[Z_] - mult*mult/l2;
                return (m->p_tr[X_]*m->p_tr[X_] + c2);
            }
            else {  //p3 is closest
                v3 temp;
                (*dot_out)*=m->ang[i*4+2];
                temp.r[0] = m->p_tr[X_] ;
                temp.r[1] = m->p_tr[Y_] - m->shP_2[Y_];
                temp.r[2] = m->p_tr[Z_] - m->shP_2[Z_];
                return LEN2(temp);
            }
        }
    }

    if(E23 < 0) {  //left side p2p3, right side of other
        v3 temp1, temp2;

        temp1.r[0] =0.0;
        temp1.r[1] =m->shP_2[Y_];
        temp1.r[2] =m->shP_2[Z_] - m->shP_1[Z_];

        temp2.r[0] =m->p_tr[X_] ;
        temp2.r[1] =m->p_tr[Y_] ;
        temp2.r[2] =m->p_tr[Z_] - m->shP_1[Z_];

        mult = temp1.r[1]*temp2.r[1] + temp1.r[2]*temp2.r[2];//DOT(temp1, temp2);
        if(mult < 0) {  //p2 is closest
            (*dot_out)*=m->ang[i*4+1];
            return LEN2(temp2);
        }
        else {
            double l2= temp1.r[1]*temp1.r[1] + temp1.r[2]*temp1.r[2];//LEN2(temp1);
            if (mult <= l2) {     //p2p3 is closest
                double c2 = temp2.r[1]*temp2.r[1] + temp2.r[2]*temp2.r[2] - mult*mult/l2;
                return (m->p_tr[X_]*m->p_tr[X_] + c2);
            }
            else {  //p3 is closest
                v3 temp;
                (*dot_out)*=m->ang[i*4+2];
                temp.r[0] = m->p_tr[X_] ;
                temp.r[1] = m->p_tr[Y_] - m->shP_2[Y_];
                temp.r[2] = m->p_tr[Z_] - m->shP_2[Z_];
                return LEN2(temp);
            }
        }
    }
}

double  modelDistP_SSE(model* m, v3* point)
{


     double b[] = {point->r[0],point->r[1],point->r[2],1};
    for (int i=0; i<m->n_m_tris;i++)
    {
        dot4x4(&(m->matr_1[i*16]), b, &(m->p_tr[i*4])); //p_tr is now xyzw in local CS
    }


    double min_dist =1e10;
    double sum_dot=0.0;//(0,0,0);
    double delta=1e-5*1e-5;
    for (int i=0; i<m->n_m_tris;i++)
    {

        double  dot_n;
        double l= dist_sse(m,i,&dot_n);//distP_tri_optimized(&(m->m_tris[i]),point,&dot_n);

        if (min_dist>=l-delta)
        {
            if  (fabs(l-min_dist)>=delta)
                sum_dot=dot_n;
            else
                sum_dot+=dot_n;

            min_dist=l;
        }
    }
    return sqrt(min_dist)*(2.0*(sum_dot>0.0)-1.0);


}
#endif // CALC_DISTANCE_H
