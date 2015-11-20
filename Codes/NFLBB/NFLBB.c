/*
Parameter Sweep for Negative feed-back loop circuit with buffering node.
*/


#include<stdio.h>
#include<math.h>
void set(float a[],int i,int flag)  // to initialize values in loops in makefile
{
  if (flag)
  {
   a[i]=.001;
  }
  else
  {
      a[i]=.1;
  }
}
void makefile()   //creates the parameterset file
{
  float k[6]={.1,.1,.1,.1,.1,.1},m[6]={.001,.001,.001,.001,.001,.001};
  int i,j,p=0;
  float x;
  FILE *f;
  f=fopen("parameterc.txt","w");
  long long int count1=0,count2=0,count3=0,count4=0;
  for(;k[0]<11;k[0]*=10)
  {
    set(k,1,0);
    for(;k[1]<11;k[1]*=50)
    {
      set(k,2,0);
      for(;k[2]<11;k[2]*=50)
      {
	set(k,3,0);
	for(;k[3]<11;k[3]*=50)
	{
	  set(k,4,0);
	  for(;k[4]<11;k[4]*=50)
	  {
	    set(k,5,0);
	    for(;k[5]<11;k[5]*=50)
	    {
	     set(m,0,1);
	     for(;m[0]<101;m[0]*=100)
	     {
	       set(m,1,1);
	       for(;m[1]<101;m[1]*=100)
	       {
		 set(m,2,1);
		 for(;m[2]<101;m[2]*=100)
		 {
		   set(m,3,1);
		   for(;m[3]<101;m[3]*=100)
		   {
		     set(m,4,1);
		     for(;m[4]<101;m[4]*=100)
		     {
		       set(m,5,1);
		       for(;m[5]<101;m[5]*=100)
		       {
			 //write into file
			 for(i=0;i<6;i++)
			 {
			   fprintf(f,"%f ",k[i]);
			 }
			 for(j=0;j<6;j++)
			 {
			   fprintf(f,"%f ",m[j]);
			}
			fprintf(f,"\n");
			if(p==0)
			{count1++;p++;}
			else if(p==1)
			{count2++;p++;}
			else if(p==2)
			{
			  count3++;p++;
			}
			else
			{
			  count4++;p++;p=p%4;
			}
			
		      }
		    }
		  }
		}
	      }
	     }
	    }
	    
	  }
    
	}
    
    }
    
    }
  }
  
  fclose(f);
  printf("%lld  %lld %lld %lld",count1,count2,count3,count4);
}




void readtest()
{
  FILE *f;
  f=fopen("parameterc.txt","r");
  int i;
  float x;
  for(i=0;i<24;i++)
  {
    fscanf(f,"%f",&x);
    printf("%f \n",x+1);
  }
  
}

void solver(float k[6],float K[6],FILE *f) //takes a parameter set and solves the system
{
	//printf("\nCall to solver");
	double a=0.5,b=.5,c=.5,I=.5,Fa=.5,Fb=.5,da=0,db=0,dc=0,max=0;
  	double dt=0.001;
	long double stop;
	int j;
	stop = pow(10,-11);
	
	do		//first equilibrium
	{
	  
	  da=dt*( I*k[0]*(1-a)/(1-a+K[0]) - Fa*k[1]*a/(a+K[1]));
	  db=dt*( c*k[2]*(1-b)/(1-b+K[2]) - Fb*k[3]*b/(b+K[3]));
	  dc=dt*( a*k[4]*(1-c)/(1-c+K[4]) - b*k[5]*c/(c+K[5]));
	  a=a+da;
	  b=b+db;
	  c=c+dc;
	}while(da >stop || db> stop || dc>stop);
	
	
	I = 0.6;
	
	do		//second equilibrium
	{
	  
	  da=dt*( I*k[0]*(1-a)/(1-a+K[0]) - Fa*k[1]*a/(a+K[1]));
	  db=dt*( c*k[2]*(1-b)/(1-b+K[2]) - Fb*k[3]*b/(b+K[3]));
	  dc=dt*( a*k[4]*(1-c)/(1-c+K[4]) - b*k[5]*c/(c+K[5]));
	  a=a+da;
	  b=b+db;
	  c=c+dc;
	  if(c > max)
	    max=c;
	  
	}while(da > stop || db> stop || dc>stop);
	
	
	
	
	double sens,pres;
	sens= fabs(((max/.5)-1)/.2);
	pres= fabs(.2/(c/.5-1));
		
	if(sens>1&&pres>10) //add the parameter set to file
	{
	 
	  for(j=0;j<6;j++)
          {
	   fprintf(f,"%f ",k[j]);
	  }
	  for(j=0;j<6;j++)
	  {
	   fprintf(f,"%f ",K[j]);
	  }
	  fprintf(f,"\n");
			 
	}
	

		
	
}


void main()
{
  //makefile();
  FILE *f,*f1;
  f=fopen("parameterc.txt","r");
  f1=fopen("validset.txt","w");
  int i;
  float k[6],m[6];
  
  while(feof(f)==0)//read parameterset and solve
    {
      //printf("\nIn loopin main");
      for(i=0;i<6;i++)
      {
	fscanf(f,"%f",&k[i]);
      }
      for(i=0;i<6;i++)
      {
	fscanf(f,"%f",&m[i]);
      }
    
      solver(k,m,f1);
      
    }//end of loop
  
  //readtest();
  
  fclose(f);
  fclose(f1);
  
}


