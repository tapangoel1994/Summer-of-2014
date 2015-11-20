#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define NODES 3
#define EDGES 9
#define TYPES_OF_EDGES 3

void initialize(int a[],int n) //set initial value to circuit string
{
  int i=0;
  while(i<n)
  {
    a[i++]=-1;
    
  }
  
  return ;
}

void print(int a[],int n)  //print the string
{
  int i;
  for(i=0;i<n;i++)
    printf("%d  ",a[i]);
  printf("\n");
  
  
 return ;
}




void perm2 (int arr[],int edges,int nodes)
{
	
	int flag=1,i=0,min=-1,max=2,j,sum,index;
		
	
	
	while(flag==1)
	{
	
		arr[i]++;
		if(arr[i]==max)
		{
			arr[i]=min;
			i++;
		}
		else
		{
			flag=0;
		}
		
	}
	
	index=edges;
	
	while(index<edges+nodes)
	{
		sum=0;
		j=nodes*(index-edges);
		while(j<nodes*(index-edges)+nodes)
		{
			sum+=arr[j++];
		}
		
		if(sum>0)
			arr[index]=-1;
		else if(sum<0)
			arr[index]=1;
		else
			arr[index]=0;
		
		index++;
	}
	
	//print(arr,edges+nodes);

}






int main()
{
	long int i,edges,nodes,delay;
	FILE *f;
	int j,no_of_edges,p,no_of_sets;
	char str[10];
	edges=EDGES;
	nodes=NODES;
	double constants,temp;
	int a[edges+nodes];
	
	
	initialize(a,edges+nodes);
	srand((unsigned int)time(NULL));
	
	for(i=0;i<pow(TYPES_OF_EDGES,EDGES);i++)
	{
		no_of_edges=0;
		no_of_sets=10;
		
		perm2(a,edges,nodes);
	
		sprintf(str,"%ld",i);
		f=fopen(str,"w");
	
		for(j=0;j<nodes+edges;j++)
		{
			fprintf(f,"%d ",a[j]);
			if(a[j]!=0)
				no_of_edges++;
		}
		
		fprintf(f,"\n");
		
		while(no_of_sets)
		{
	
			for(p=0;p<no_of_edges*2+2;p++)
			{
								
				
				//constants=log unif(.001 - 100) if p%2==1 and constants=log unif(.1 - 10) if p%2==0	
				if(p%2==0)
				{
					
					temp=(double)rand()/RAND_MAX;
					temp=temp*3 - 1 ;
					constants=pow(10,temp);
									
				}
				
				
				else
				{
					temp=(double)rand()/RAND_MAX;
					temp=temp*5 - 3; 
					constants=pow(10,temp);
					
				}
					
				fprintf(f,"%lf ",constants);
			
			}
			fprintf(f,"\n");
			
			no_of_sets--;
		
		}
		
		fclose(f);
		
	}

  
  	return 0;
}