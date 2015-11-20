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


double *x,*xnew,*k2,*k3,*k4,*k5,*k6,I;
int nodes=3,edges=9,*circ;
double *constants;

void eval1(double t, double *x_init, double *f)
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
			else if(circ[i*nodes+j]==-1)
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
	f[0]+=constants[p++]*I*(1-x_init[0])/(constants[p++]+1-x_init[0]);
	return ;
}




double rkstep(double t,double *x,double *xdot,double h,double *xout)
{
// 	printf("3\n");
    int j;
    double xtemp[nodes],xerr,errmax=0;
    for (j=0;j<nodes;j++) xtemp[j]=x[j]+b21*h*xdot[j];
    eval1(t,xtemp,k2);
// 	printf("3.1\n");
    for (j=0;j<nodes;j++) xtemp[j]=x[j]+h*(b31*xdot[j]+b32*k2[j]);
    eval1(t,xtemp,k3);
// 	printf("3.2\n");
    for (j=0;j<nodes;j++) xtemp[j]=x[j]+h*(b41*xdot[j]+b42*k2[j]+b43*k3[j]);
    eval1(t,xtemp,k4);
// 	printf("3.3\n");
    for (j=0;j<nodes;j++) xtemp[j]=x[j]+h*(b51*xdot[j]+b52*k2[j]+b53*k3[j]+b54*k4[j]);
    eval1(t,xtemp,k5);
// 	printf("3.4\n");
    for (j=0;j<nodes;j++) xtemp[j]=x[j]+h*(b61*xdot[j]+b62*k2[j]+b63*k3[j]+b64*k4[j]+b65*k5[j]);
    eval1(t,xtemp,k6);
// 	printf("3.5\n");
    for (j=0;j<nodes;j++)
        {
        xout[j]=x[j]+h*(c1*xdot[j]+c3*k3[j]+c4*k4[j]+c6*k6[j]);
        xerr=fabs(h*(dc1*xdot[j]+dc3*k3[j]+dc4*k4[j]+dc5*k5[j]+dc6*k6[j])*EPSINV/(fabs(x[j])+h*fabs(xdot[j])+TINY));
        if (xerr>errmax) errmax=xerr;
        }
// 	printf("4\n");
    return errmax;
}

    
void deterministicrun(double ti,double tf,double *x0,int storeresult)
{
    int i,j,nsteps,run;
    double length=0,h,t,err,errx,errxmax=1;
    int flag=0,tspanflag=0;
	double xdot1[nodes];
	double c_eq=0,c_max=0,c_final=0;
	
	
    for (i=0;i<nodes;i++) x[i]=x0[i];
    h=0.001;nsteps=0;
    t=ti;
    int flageq=1;
    I=0.003;
    while (/*t<tf &&*/ flageq)
    {   
		//printf("0\n");
        eval1(t,x,xdot1);
// 		printf("1\n");
        do 
        {
            err=rkstep(t,x,xdot1,h,xnew);
// 			printf("2");
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
            if(fabs(xdot1[i])>EPSEQ)
                flageq++;
        }
        printf("%lf  ",t);
        for (i=0;i<nodes;i++)
            printf("%lf  ",x[i]);
        for (i=0;i<nodes;i++)
            printf("%lf  ",xdot1[i]);
        printf("%lf  ",I);
        printf("\n");
    }
    
    c_eq=x[nodes-1];  //set initial equilibrium value of c
    
    printf("%lf  ",t);
    for (i=0;i<nodes;i++)
       printf("%lf  ",x[i]);
    for (i=0;i<nodes;i++)
       printf("%lf  ",xdot1[i]);
	printf("%lf  ",I);
    printf("\n");
  
  
    I=1.2*I;
    flageq=1;
    h=0.001;
    while (/*t<2*tf  &&*/ flageq)
    {   
		if(c_max<x[nodes-1])   //set cmax to max value
			c_max=x[nodes-1];
        //printf("0\n");
        eval1(t,x,xdot1);
// 		printf("1\n");
        do 
        {
            err=rkstep(t,x,xdot1,h,xnew);
// 			printf("2");
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
            if(fabs(xdot1[i])>EPSEQ)
                flageq++;
        }
        printf("%lf  ",t);
        for (i=0;i<nodes;i++)
            printf("%lf  ",x[i]);
        for (i=0;i<nodes;i++)
            printf("%lf  ",xdot1[i]);
        printf("%lf  ",I);
        printf("\n");
    }
    
    c_final=x[nodes-1]; //set value of cfinal
	
    return ;
}



int main()
{
	FILE *fp;
	long int i,j,k;
	char str[10];
	
	nodes=3;
	edges=9;
	circ=new int [nodes + edges];
	int no_of_edges;
	int no_of_constants;
	x=new double [nodes];
	k2=new double [nodes];
	k3=new double [nodes];
	k4=new double [nodes];
	k5=new double [nodes];
	k6=new double [nodes];
 	xnew=new double [nodes];
	
	
	
	for(i=10;i<11;i++)//loop to open a file in read mode and then an internal loop to run each parameter set based on circuit.
	{
		sprintf(str,"%ld",i);
		fp=fopen(str,"r");
		
		no_of_edges=0;
		
		for(j=0;j<nodes+edges;j++) // to get the circuit
		{
			fscanf(fp,"%d",&circ[j]);
			if(circ[j]!=0)
				no_of_edges++;
		}
		
		//print(circ,nodes+edges);
		
		no_of_constants=no_of_edges*2+2;
		
		constants=new double [no_of_constants];
		
		for(k=0;k<1;k++)  //iterate over  parameter sets
		{
			for(j=0;j<no_of_constants;j++)
			{
				fscanf(fp,"%lf",&constants[j]);
			}
			
			for(j=0;j<nodes;j++)
				{
					x[j]=.5;
// 					printf("%lf ",x[j]);
				}
		//	printf("\n");
			deterministicrun(0,10000,x,0);
		}
		fflush(fp);
		
// 		fclose(fp);
	
	}
	
	exit(0);
	return 0;
	
}