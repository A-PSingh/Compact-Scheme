//Question1
// compact differencing scheme for first order drivative
#include<stdio.h>
#include<math.h>
int main(void)
{
	int i,j,k,M=101,N=101,order;
	double phi[M][N],dphibydx[M][N],dphibydxinitial[M][N],P[M][N],Q[M][N],d[M][N];
	double dx,dy,xi,yi,x,y,beta,alpha,a,b,c,a1,ai,bi,b1,ci,t1,t2,t3,t4,t5,t6,temp;
	dx=10.0/(M-1); 
	dy=10.0/(N-1); 
	 FILE *fp;
	for(i=0;i<M;i++)
  	{ for(j=0;j<N;j++)
	   {
		xi=(i)*dx; yi=(j)*dy;
	phi[i][j] =sin(xi) *cos(yi);
    } }
   	printf("Enter the order of accuracy:4,6 or 8  ");
    scanf("%d",&order);
    
// One sided fourth order explicit scheme for boundary nodes and nodes adjacent to boundary nodes
   for(i=0;i<=1;i++){
     for(j=0;j<N;j++)
        {//Expression for 4th order forward differencing:
           dphibydx[i][j]= (-1*25*phi[i][j] +48*phi[i+1][j] -36*phi[i+2][j] +16*phi[i+3][j] -3*phi[i+4][j]) /(12*dx);
		}}
   for(i=M-1;i>M-3;i--)
    { for(j=0;j<N;j++)
        {//Expression for 4th backward differencing:
           dphibydx[i][j]= (25*phi[i][j] -48*phi[i-1][j] +36*phi[i-2][j] -16*phi[i-3][j] +3*phi[i-4][j]) /(12*dx);
	 	}} 
	
	if(order==4){
	printf("for fourth order accuracy\n"); 	
    alpha=0.3;
    //for fourth order of accuracy
    beta=0; a=(2.0/3)*(alpha+2); b=(1.0/3)*(4*alpha-1); c=0;	                                                    
 //TDMA
     for(i=2;i<=M-3;i++){
        if(i==2)
           for(j=0;j<N;j++){
       d[i][j]= a*(phi[i+1][j]-phi[i-1][j])/(2*dx)   +b*(phi[i+2][j]-phi[i-2][j])/(4*dx)  -alpha*dphibydx[1][j];
			    }
    	else if(i==M-3) 
           for(j=0;j<N;j++){
       d[i][j]= a*(phi[i+1][j]-phi[i-1][j])/(2*dx)   +b*(phi[i+2][j]-phi[i-2][j])/(4*dx) -alpha*dphibydx[M-2][j];
           
		       }
     	else
           for(j=0;j<N;j++){
       d[i][j]= a*(phi[i+1][j]-phi[i-1][j])/(2*dx)   +b*(phi[i+2][j]-phi[i-2][j])/(4*dx);
           }
             }
     //By writing equations in standard form: ai *dphibydx[i]=bi *dphibydx[i+1] +ci *dphibydx[i-1] +di; we get all ai, bi and ci are same
	 ai=1; bi=-1*alpha; ci=-1*alpha; a1=1; b1=-1*alpha; 
	   
  for(j=0;j<N;j++){
      P[2][j]=b1/a1; 
	  Q[2][j]=d[2][j]/a1;
   }
  for(i=3;i<=M-3;i++){
  	for(j=0;j<N;j++){
     P[i][j]=bi/(ai-ci*P[i-1][j]);
     Q[i][j]=(d[i][j]+ci*Q[i-1][j]) / (ai-ci*P[i-1][j]);	
    }}
    
	for(j=0;j<N;j++)
      dphibydx[M-3][j] =Q[M-3][j];

	for(j=0;j<N;j++) 
	{  
    for(i=M-4;i>=2;i--)
        { dphibydx[i][j] = P[i][j] *dphibydx[i+1][j] + Q[i][j];
      }	                                                                                               
   }	
  
   fp=fopen("4thorder.dat","w");
   for(i=0;i<M;i++){
     for(j=0;j<N;j++){
     	x=(i)*dx; 
     	y=(j)*dy;	
   printf("%lf	%lf	%lf\n",x,y,dphibydx[i][j]);
   fprintf(fp,"%lf	%lf	%lf\n",x,y,dphibydx[i][j]);
          }} 
          fclose(fp);
    fp=fopen("4dphibydxmidvertical.dat","w");
     for(j=0;j<N;j++){
     	y=(j)*dy;	 
   fprintf(fp,"%lf	%lf\n",dphibydx[50][j],y);
          }
          fclose(fp);
	 fp=fopen("4dphibydxmidhorizontal.dat","w");
     for(i=0;i<N;i++){
     	x=(i)*dx; 
   fprintf(fp,"%lf	%lf\n",x,dphibydx[i][50]);
          }
          fclose(fp);
    t1=dphibydx[30][50]; t2=dphibydx[60][30];  
	    
   }
   
   if(order==6){    
    //6th order
   	printf("for sixth order accuracy\n");
    
	// One sided fourth order explicit scheme for boundary nodes and nodes adjacent to boundary nodes
   for(i=0;i<=M-1;i++){
   	     if(i==0||i==1||i==2){
              for(j=0;j<N;j++)
                   {//Expression for 4th order forward differencing:
           dphibydx[i][j]= (-1*25*phi[i][j] +48*phi[i+1][j] -36*phi[i+2][j] +16*phi[i+3][j] -3*phi[i+4][j]) /(12*dx);
       // printf("%lf\n",dphibydx[i][j]);
	             	}}
         else if(i==M-1||i==M-2||i==M-3)
                { for(j=0;j<N;j++)
                       {//Expression for 4th backward differencing:
           dphibydx[i][j]= (25*phi[i][j] -48*phi[i-1][j] +36*phi[i-2][j] -16*phi[i-3][j] +3*phi[i-4][j]) /(12*dx);
       // printf("%lf\n",dphibydx[i][j]);
	               	}} 
        	else 
	           dphibydx[i][j]= 0.0;	 
		          }
		 
	  for(i=0;i<=M-1;i++){
	  	      for(j=0;j<N;j++){
	             dphibydxinitial[i][j]=dphibydx[i][j];
	             //if(i==2)
				 //printf("%lf\n",dphibydx[i][j]);
                       }}
    alpha=0.3;	 
    //for sixth order of accuracy
    beta=0; a=(1.0/6)*(alpha+9); b=(1.0/15)*(32*alpha-9); c=(1.0/10)*(-3*alpha +1);
                                                                                                                                                                        
//Gauss seidel method
 float error=1.0,temp;
 
 do{  error=0.0;
  for (i=3;i<=M-4;i++){
           for(j=0;j<N;j++){
                   
  	dphibydx[i][j] =-alpha*dphibydx[i+1][j] -beta*dphibydx[i+2][j]-beta*dphibydx[i-2][j]-alpha*dphibydx[i-1][j] 
	             +a*(phi[i+1][j]-phi[i-1][j])/(2*dx) +b*(phi[i+2][j]-phi[i-2][j])/(4*dx);
       
           error=error+fabs(dphibydxinitial[i][j]-dphibydx[i][j]);
           dphibydxinitial[i][j]=dphibydx[i][j];
                   }} 
 }while(error>1e-8);     
   fp=fopen("6gauss.dat","w");                                                                                        
    for (i=0;i<M;i++){
       for(j=0;j<N;j++){
	   x=(i)*dx; 
     	y=(j)*dy;	
   printf("%lf	%lf	%lf\n",x,y,dphibydx[i][j]);
   fprintf(fp,"%lf	%lf	%lf\n",x,y,dphibydx[i][j]);
        }}
               fclose(fp); 
			   
	fp=fopen("6dphibydxmidvertical.dat","w");
     for(j=0;j<N;j++){
     	y=(j)*dy;	 
   fprintf(fp,"%lf	%lf\n",dphibydx[50][j],y);
          }
          fclose(fp);
	 fp=fopen("6dphibydxmidhorizontal.dat","w");
     for(i=0;i<N;i++){
     	x=(i)*dx; 
   fprintf(fp,"%lf	%lf\n",x,dphibydx[i][50]);
          }
          fclose(fp);		    
	t3=dphibydx[30][50]; t4=dphibydx[60][30];
		                                                                                
  }
 if(order==8){
                  
	//8th order
	printf("for eighth order accuracy\n");	   
	// One sided fourth order explicit scheme for boundary nodes and nodes adjacent to boundary nodes
   for(i=0;i<=M-1;i++){
   	     if(i==0||i==1||i==2){
              for(j=0;j<N;j++)
                   {//Expression for 4th order forward differencing:
           dphibydx[i][j]= (-1*25*phi[i][j] +48*phi[i+1][j] -36*phi[i+2][j] +16*phi[i+3][j] -3*phi[i+4][j]) /(12*dx);
       // printf("%lf\n",dphibydx[i][j]);
	             	}}
         else if(i==M-1||i==M-2||i==M-3)
                { for(j=0;j<N;j++)
                       {//Expression for 4th backward differencing:
           dphibydx[i][j]= (25*phi[i][j] -48*phi[i-1][j] +36*phi[i-2][j] -16*phi[i-3][j] +3*phi[i-4][j]) /(12*dx);
       // printf("%lf\n",dphibydx[i][j]);
	               	}} 
        	else 
	           dphibydx[i][j]= 0.0;	 
		          }
		 
	  for(i=0;i<=M-1;i++){
	  	      for(j=0;j<N;j++){
	             dphibydxinitial[i][j]=dphibydx[i][j];
	             //if(i==2)
				 //printf("%lf\n",dphibydx[i][j]);
                       }}
    alpha=0.3;
    beta=(1.0/20)*(-3+8*alpha); a=(1.0/6)*(12-7*alpha); b=(1.0/150)*(568*alpha-183); c=(1.0/50)*(9*alpha -4);   
// General five point formulation for the approximation of first order derivative
                                                                                                                                                                                            
//Gauss seidel method
  double error=1.0 ;
 
 do{  error=0.0;
  for (i=3;i<=M-4;i++){
           for(j=0;j<N;j++){
                   
  	dphibydx[i][j] =-alpha*dphibydx[i+1][j] -beta*dphibydx[i+2][j]-beta*dphibydx[i-2][j]-alpha*dphibydx[i-1][j] 
	             +a*(phi[i+1][j]-phi[i-1][j])/(2*dx) +b*(phi[i+2][j]-phi[i-2][j])/(4*dx);
       
           error=error+fabs(dphibydxinitial[i][j]-dphibydx[i][j]);
           dphibydxinitial[i][j]=dphibydx[i][j];
                   }}
 }while(error>1e-8);  
   fp=fopen("8gauss.dat","w");                                                                                        
    for (i=0;i<M;i++){
       for(j=0;j<N;j++){
	   x=(i)*dx; 
     	y=(j)*dy;	
   printf("%lf	%lf	%lf\n",x,y,dphibydx[i][j]);
   fprintf(fp,"%lf	%lf	%lf\n",x,y,dphibydx[i][j]);
        }}
               fclose(fp);
               
    fp=fopen("8dphibydxmidvertical.dat","w");
     for(j=0;j<N;j++){
     	y=(j)*dy;	 
   fprintf(fp,"%lf	%lf\n",dphibydx[50][j],y);
          }
          fclose(fp);
	 fp=fopen("8dphibydxmidhorizontal.dat","w");
     for(i=0;i<N;i++){
     	x=(i)*dx; 
   fprintf(fp,"%lf	%lf\n",x,dphibydx[i][50]);
          }
          fclose(fp);           
	t5=dphibydx[30][50]; t6=dphibydx[60][30]; 

}
	//ANALYTICAL METHOD
			     
	double Analytical[M] [N];
	fp=fopen("Analytical.dat","w"); 
	for(i=0;i<M;i++)
  	{ for(j=0;j<N;j++)
	   {
		xi=(i)*dx;
   	   yi=(j)*dy;
    	Analytical[i] [j]=cos(xi) *cos(yi);
        printf("%lf	%lf	%lf\n",xi,yi,Analytical[i][j]);
        fprintf(fp,"%lf	%lf	%lf\n",xi,yi,Analytical[i][j]);
    }}		                                                                                
    fclose(fp);  
	fp=fopen("Analyticalmidvertical.dat","w");
     for(j=0;j<N;j++){
     	yi=(j)*dy;	 
   fprintf(fp,"%lf	%lf\n",Analytical[50][j],yi);
          }
          fclose(fp);
	 fp=fopen("Analyticalmidhorizontal.dat","w");
     for(i=0;i<N;i++){
     	xi=(i)*dx; 
   fprintf(fp,"%lf	%lf\n",xi,Analytical[i][50]);
          }
          fclose(fp); 
    if(order==4)
	{
	printf("At (3,5)\t(6,3)\n");	  	        	    
	printf("%lf\t%lf\n",t1,t2);      
	printf("By Analytical method\n");	    
	printf("%lf\t%lf\n",Analytical[30][50],Analytical[60][30]);	 
    }
    return(0);
}
 
