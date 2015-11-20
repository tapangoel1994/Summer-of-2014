#include<stdio.h>

int main()
{
    float a=.5,as=.5,f=1,af=0,aas=0,I=0,aI=0,an=0,asn=0,fn=0,afn=0,aasn=0,aIn=0;
    float k1=.01,k_1=.1,k2=1,k_2=.1,k3=1,k_3=.01,k4=.1,k_4=.01,k5=.1,k_5=0,k6=1,k_6=.1,k7=.1,k_7=.01,k8=0.1;
    float dt=0.001;
    int i;
    float diff=0;
    float b[20000][7];
    for (i=0;i<20000;i++)
    {


        an=a+dt*(-k1*a*as + k_1*aas +2*k3*aas - k_3*a*a -k4*a*f +k_4*af +k7*aI -k_7*a*I);
        asn=as+dt*(-k1*a*as + k_1*aas + 2*k2*aas -k_2*as*as +k5*af -k_5*as*f - k6*as*I +k_6*aI);
        aasn=aas+dt*(k1*a*as -k_1*aas +0.5*k_2*as*as -k2*aas +0.5*k_3*a*a - k3*aas);
        afn = af+dt*(k4*a*f -k_4*af -k5*af +k_5*as*f);
        fn= f + (af-afn) +dt*k8*(1-f);
        aIn= aI + dt*(k6*as*I -k_6*aI - k7*aI + k_7*a*I);

        a=an;as=asn;aas=aasn;af=afn;f=fn;aI=aIn;
        b[i][0]=a;
        b[i][1]=as;
        b[i][2]=aas;
        b[i][3]=af;
        b[i][4]=f;
        b[i][5]=aI;
        b[i][6]=a+as+2*aas+af+aI;


    }

    FILE *fp;
    int j;
    fp=fopen("data.txt","w");
    for(i=0;i<20000;i++)
    {
        for ( j=0;j<7;j++)
        {
            putc(b[i][j],fp);
        }
    }


}
