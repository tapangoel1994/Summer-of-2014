#include<stdio.h>
#include<stdlib.h>
#include<math.h>

using namespace std;

// various constants required by integrator
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4   
#define TINY 1.0e-30
#define MAXSTP 10000
#define EPSINV 1.0e11       //modified from 1.0e5
#define EPSEQ 1e-5
#define EXT .5

const double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,
    b21=0.2,b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42=-0.9,b43=1.2,
    b51=-11.0/54.0,b52=2.5,b53=-70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,
    c1=37.0/378.0,c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,dc4=c4-13525.0/55296.0,
    dc5=-277.0/14336.0,dc6=c6-0.25;

// ---------------------------


int n,nedge,l;
double *x,*xdot,*xnew,*k2,*k3,*k4,*k5,*k6,I;
double *con1,*con2,*conm1,*conm2;
int *from,*to;

void eval1(double t, double *x_init, double *f, int *circ,float *constants )
{
	int j;
	int i;
	int p=0;
	for(j=0;j<nodes;j++)
	{
		f[j]=0;
		for(i=0;i<nodes;i++)
		{
			if(circ[i*nodes+j]==-1)
			{
				f[j]-= constants[p++]*x_init[i]*x_init[j]/(constants[p++]+x_init[j]);
			}
			else if(circ[i*n+j]==-1)
			{
				f[j]+= constants[p++]*x_init[i]*(1-x_init[j])/(constants[p++]+1-x_init[j]);
			}
		}
	}
    for(j=nodes*nodes;j<nodes*nodes+nodes;j++)
	{
		if(circ[j]==-1)
		{
			f[j]-= constants[p++]*EXT*x_init[j]/(constants[p++]+x_init[j]);   
		}
		else if(circ[j]==1)
		{
			f[j]+=constants[p++]*EXT*(1-x_init[j])/(constants[p++]+1-x_init[j]);   
		}
	}
}



void eval(double t,double *xless,double *f) // new equations with cost for RM
{
    int j;
	
    double dB, dBs;
//     f[0]=5*0.5*xless[0]/(0.001+xless[0]);
//     f[l]=-5*0.5*xless[0]/(0.001+xless[0]);
//     f[1]=5*0.5*xless[1]/(0.001+xless[1]);
//     f[l+1]=-5*0.5*xless[1]/(0.001+xless[1]);
    for(j=0;j<l;j++)
    {
        f[j]=-0.1*xless[j] + 0.1*xless[l+j];
        f[l+j]=0.1*xless[j] - 0.1*xless[l+j];
    }
    for(j=0;j<nedge;j++)
    {
        dBs= xless[from[j]]*(-1*con1[j]*conm2[j]*xless[l+to[j]]+ con2[j]*conm1[j]*xless[to[j]])/(conm1[j]+conm2[j])*(con1[j]*xless[l+to[j]]+con2[j]*xless[to[j]]+conm2[j]+conm1[j])/(con1[j]*xless[l+to[j]]+con2[j]*xless[to[j]]+conm2[j]+conm1[j]+con1[j]*xless[from[j]] );
        
        dB= xless[from[j]]*(con1[j]*conm2[j]*xless[l+to[j]]-con2[j]*conm1[j]*xless[to[j]])/(conm1[j]+conm2[j])*(con1[j]*xless[l+to[j]]+con2[j]*xless[to[j]]+conm2[j]+conm1[j])/(con1[j]*xless[l+to[j]]+con2[j]*xless[to[j]]+conm2[j]+conm1[j]+con2[j]*xless[from[j]] );
        
        f[to[j]]+= dB;
        f[to[j]+l]+= dBs;
        
        f[from[j]]+= -(xless[from[j]]*con1[j])/(con1[j]*xless[to[j]+l]+con2[j]*xless[to[j]]+conm2[j]+conm1[j])*dBs  -(xless[from[j]]*con2[j])/(con1[j]*xless[to[j]+l]+con2[j]*xless[to[j]]+conm2[j]+conm1[j])*dB;
    }
    f[0]+=10*I*xless[l]/(0.01+xless[l]);
    f[l]+=-10*I*xless[l]/(0.01+xless[l]);
}

double rkstep(double t,double *x,double *xdot,double h,double *xout)
{
    int j;
    double xtemp[nodes],xerr,errmax=0;
    for (j=0;j<nodes;j++) xtemp[j]=x[j]+b21*h*xdot[j];
    eval(t,xtemp,k2);
    for (j=0;j<nodes;j++) xtemp[j]=x[j]+h*(b31*xdot[j]+b32*k2[j]);
    eval(t,xtemp,k3);
    for (j=0;j<nodes;j++) xtemp[j]=x[j]+h*(b41*xdot[j]+b42*k2[j]+b43*k3[j]);
    eval(t,xtemp,k4);
    for (j=0;j<nodes;j++) xtemp[j]=x[j]+h*(b51*xdot[j]+b52*k2[j]+b53*k3[j]+b54*k4[j]);
    eval(t,xtemp,k5);
    for (j=0;j<nodes;j++) xtemp[j]=x[j]+h*(b61*xdot[j]+b62*k2[j]+b63*k3[j]+b64*k4[j]+b65*k5[j]);
    eval(t,xtemp,k6);
    for (j=0;j<nodes;j++)
        {
        xout[j]=x[j]+h*(c1*xdot[j]+c3*k3[j]+c4*k4[j]+c6*k6[j]);
        xerr=fabs(h*(dc1*xdot[j]+dc3*k3[j]+dc4*k4[j]+dc5*k5[j]+dc6*k6[j])*EPSINV/(fabs(x[j])+h*fabs(xdot[j])+TINY));
        if (xerr>errmax) errmax=xerr;
        }
    return errmax;
}

    
void deterministicrun(double ti,double tf,double *x0,int storeresult)
{
    int i,j,nsteps,run;
    double length=0,h,t,err,errx,errxmax=1;
    int flag=0,tspanflag=0;

    for (i=0;i<nodes;i++) x[i]=x0[i];
    h=0.001;nsteps=0;
    t=ti;
    int flageq=1;
    I=0.003;
    while (/*t<tf &&*/ flageq)
    {   
        eval(t,x,xdot);
        do 
        {
            err=rkstep(t,x,xdot,h,xnew);
            //printf("Err %.20lf  ",err);
            if (err>1.0/* && h>0.00028*/)   //pritish change err>1.0 to err>0.1
            {
                double htemp=SAFETY*h*pow(err,PSHRNK);
                if (htemp<0.1*h) htemp=0.1*h;
                h=htemp;tspanflag=0;
            }
            else 
            {
                t+=h;//nsteps++;
//                 if (storeresult==1) {for (i=0;i<n;i++) {printf("%lf ",xnew[i]);}printf("\n");} 
                for (i=0;i<nodes;i++)
                        {
                    x[i]=xnew[i];
                    //if (x[i]<0) x[i]=0;
                    }
                if (err>ERRCON) h=SAFETY*h*pow(err,PGROW);
                else h*=5.0;
            }
        }while (err>1.0);    //pritish change err>1.0 to err>0.1
        
        flageq=0;
        for (i=0;i<nodes;i++)
        {
            if(fabs(xdot[i])>EPSEQ)
                flageq++;
        }
        printf("%lf\t",t);
        for (i=0;i<nodes;i++)
            printf("%lf\t",x[i]);
        for(i=0;i<l;i++)
            printf("%lf\t",x[i]+x[l+i]);
        for (i=0;i<nodes;i++)
            printf("%lf\t",xdot[i]);
        printf("%lf\t",I);
        printf("\n");
    }
    
        printf("%lf\t",t);
        for (i=0;i<nodes;i++)
            printf("%lf\t",x[i]);
        for(i=0;i<l;i++)
            printf("%lf\t",x[i]+x[l+i]);
        for (i=0;i<nodes;i++)
            printf("%lf\t",xdot[i]);
        printf("%lf\t",I);
        printf("\n");
    
    
    I=1.2*I;
    flageq=1;
    h=0.001;
    while (/*t<2*tf  &&*/ flageq)
    {   
        eval(t,x,xdot);
        do 
        {
            err=rkstep(t,x,xdot,h,xnew);
            //printf("Err %.20lf  ",err);
            if (err>1.0/* && h>0.00028*/)   //pritish change err>1.0 to err>0.1
            {
                double htemp=SAFETY*h*pow(err,PSHRNK);
                if (htemp<0.1*h) htemp=0.1*h;
                h=htemp;tspanflag=0;
            }
            else 
            {
                t+=h;//nsteps++;
//                 if (storeresult==1) {for (i=0;i<n;i++) {printf("%lf ",xnew[i]);}printf("\n");} 
                for (i=0;i<nodes;i++)
                        {
                    x[i]=xnew[i];
                    //if (x[i]<0) x[i]=0;
                    }
                if (err>ERRCON) h=SAFETY*h*pow(err,PGROW);
                else h*=5.0;
            }
        }while (err>1.0);    //pritish change err>1.0 to err>0.1
        
        flageq=0;
        for (i=0;i<nodes;i++)
        {
            if(fabs(xdot[i])>EPSEQ)
                flageq++;
        }
        printf("%lf\t",t);
        for (i=0;i<nodes;i++)
            printf("%lf\t",x[i]);
//         for(i=0;i<l;i++)
//             printf("%lf\t",x[i]+x[l+i]);
//         for (i=0;i<n;i++)
//             printf("%lf\t",xdot[i]);
        printf("%lf\t",I);
        printf("\n");
    }
}


int main()
{
    int nodes=3;
    int i;
    
    for(i=0;i<nodes;i++)
    {
        x[i]=0.5;
    }

    con1= new double [nedge];
    con2= new double [nedge];
    conm1= new double [nedge];
    conm2= new double [nedge];
    
    from= new int [nedge];
    to= new int [nedge];      
    
    for(i=0;i<nedge;i++)
    {
//         printf("give one set\n");
        scanf("%d%d%lf%lf%lf%lf",&from[i],&to[i],&con1[i],&conm1[i],&con2[i],&conm2[i]);
//         printf("one set taken\n");
    }
    deterministicrun(0,10000,x,0);
    return 0;
}
