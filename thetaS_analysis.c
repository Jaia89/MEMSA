#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<float.h>

 int *k1,*kc1,*q1,**knn1,*k2,*kc2,*q2,**knn2;



void shuffle(int *array, size_t n) {    
    struct timeval tv;
    gettimeofday(&tv, NULL);
    int usec = tv.tv_usec;
    srand48(usec);


    if (n > 1) {
        size_t i;
        for (i = n - 1; i > 0; i--) {
            size_t j = (unsigned int) (drand48()*(i+1));
            int t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

/* Define degree classes for the nodes*/
int classes(int *k, int maxk, int T, int *kc){
   
	int i,j,m;
        double M;

    M=log((double)maxk)/log(2);  

    for(m=1;m<ceil(M)+1;m++){
        
    for (i=0;i<T;i++){   
    
    if((k[i]>pow((m-1),2))&&(k[i]<=pow(m,2))){
        
    kc[i]=m;}
             
       }
    }
    return(ceil(M));
}


/* Compute the entropy of the network block model*/
double entropy(int T, int **knn, int *k, int *kc, int *q, int maxk, int maxkc, int maxq, int stirling_mode){ 
       
        int i,j,l,m,c,f,h,n,ci,cl,cf,cn,d,**index,*n1_kq,**n_kq,**l_kqk1q1;      
        double fact1,fact2,fact3,S,A,B,L,c1;

        n_kq=(int**)calloc(maxkc+1,sizeof(int*));  
        
        for(i=0;i<maxkc+1;i++){
              n_kq[i]=(int*)calloc(maxq+1,sizeof(int));}
       
        index=(int**)calloc(maxkc+1,sizeof(int*));  
        
        for(i=0;i<maxkc+1;i++){
              index[i]=(int*)calloc(maxq+1,sizeof(int));}

      
      if (index==NULL){printf("error: array index not allocated!\n");}
      if (n_kq==NULL){printf("error: array n_kq not allocated!\n");}
      

      /*Block structure definition*/
 
      for (i=0;i<T;i++){
                     //printf("i=%d k[i]=%d, q[i]=%d, nkq=%d, \n",i, k[i],q[i], n_kq[k[i]][q[i]]);
                     n_kq[kc[i]][q[i]]+=1;}
    
      d=0;
       
      for (i=0;i<maxkc+1;i++){	
                 for(j=0;j<maxq+1;j++){

                 if (n_kq[i][j]>0){                 
                 d+=1;                  
                 index[i][j] = d;

                 }
             }
          }
     


           n1_kq=(int*)calloc(d+1,sizeof(int));

           if (n1_kq==NULL){printf("error: array n1_kq not allocated\n");}
        
               for (i=0;i<maxkc+1;i++){	
                 for(j=0;j<maxq+1;j++){
                 
                 if (n_kq[i][j]>0){
                 
                 n1_kq[index[i][j]]+=n_kq[i][j];
                 
               }
             }
            }



      for(i=0;i<maxkc+1;i++){free(n_kq[i]);}
      free(n_kq);


            l_kqk1q1=(int**)calloc(d+1,sizeof(int*));
       
            if (l_kqk1q1==NULL){printf("error: array l_kqk1q1 not allocated\n");}

            for(i=0;i<d+1;i++){
	      l_kqk1q1[i]=(int*)calloc(d+1,sizeof(int));

            //if (l_kqk1q1[i]==NULL){printf("error,%d  \n",i); exit(0);}
            if (l_kqk1q1[i]==NULL){printf(" Memory error! Please use degree-classes mode\n"); exit(0);}                                                           

            }




         for (i=0;i<T;i++){
                          
               if(k[i]>0){
                            
                   for(j=0;j<k[i];j++){
                         
                        if(knn[i][j]>i){                        
                        l_kqk1q1[index[kc[i]][q[i]]][index[kc[knn[i][j]]][q[knn[i][j]]]] += 1;}
              }
            }
          }

       
      for(i=0;i<maxkc+1;i++){free(index[i]);}
      free(index);
 


       // Entropy computed with Stirling approximation 
       if(stirling_mode==1){       

        
       for (l=0;l<d+1;l++){	
                     for(m=0;m<d+1;m++){
             

          if (l_kqk1q1[l][m]>0){  

          if (l==m){

                fact1=ceil((n1_kq[l]*(n1_kq[m]-1))/2);
                fact2=ceil((n1_kq[l]*(n1_kq[m]-1))/2)-l_kqk1q1[l][m];
                fact3=l_kqk1q1[l][m];

                if(ceil((n1_kq[l]*(n1_kq[m]-1))/2)<=1){
                A=0.;
                }else{
                A=(fact1*log(fact1)-fact1);
                }

                if(ceil((n1_kq[l]*(n1_kq[m]-1))/2)==2){
                A=log(2.);}
                
                
                if(ceil((n1_kq[l]*(n1_kq[m]-1))/2)-l_kqk1q1[l][m]<=1){
                B=0.;
                }else{
                B=(fact2*log(fact2)-fact2);}

                if(ceil((n1_kq[l]*(n1_kq[m]-1))/2)-l_kqk1q1[l][m]==2){
                B=log(2.);}


                if(l_kqk1q1[l][m]<=1){
                L=0.;
                }else{
                L=(fact3*log(fact3)-fact3);}

                if(l_kqk1q1[l][m]==2){
                L=log(2.);}

                
                S+=A-B-L;

                }else{

                fact1= n1_kq[l]*n1_kq[m];
                fact2= n1_kq[l]*n1_kq[m]-l_kqk1q1[l][m];
                fact3= l_kqk1q1[l][m];
               
                
                if(n1_kq[l]*n1_kq[m]<=1)
                {A=0.;}else{A=(fact1*log(fact1)-fact1);}
                if(n1_kq[l]*n1_kq[m]==2)
                {A=log(2.);}

                if(n1_kq[l]*n1_kq[m]-l_kqk1q1[l][m]<=1)
                {B=0.;}else{B=(fact2*log(fact2)-fact2);}
                if(n1_kq[l]*n1_kq[m]-l_kqk1q1[l][m]==2)
                {B=log(2.);}


                if(l_kqk1q1[l][m]<=1){
                L=0.;
                }else{
                L=(fact3*log(fact3)-fact3);}

                if(l_kqk1q1[l][m]==2){
                L=log(2.);}
                

                S+=A-B-L;

             }
          
             }else{S+=0.;}


          }
         }
    

         }else{
              

         // Entropy computed without Stirling approximation

         for (l=0;l<d+1;l++){	
                     for(m=0;m<d+1;m++){
             
          if (l_kqk1q1[l][m]>0){  

          if (l==m){

                fact1=ceil((n1_kq[l]*(n1_kq[m]-1))/2);
                fact2=ceil((n1_kq[l]*(n1_kq[m]-1))/2)-l_kqk1q1[l][m];
                fact3=l_kqk1q1[l][m];

                if(ceil((n1_kq[l]*(n1_kq[m]-1))/2)-l_kqk1q1[l][m]<=1){
 
                A=0.;
                }else{
                c1=1;            
                for(c = fact2; c <= fact1; c++)
                 {c1 = c1 * c;}
                A=log(c1);
                }
  

                if(l_kqk1q1[l][m]<=1){
                L=0.;
                }else{
                c1=1;            
                for(c = 1; c <= fact3; c++)
                 {c1 = c1 * c;}
                L=log(c1);
                }
 
                
                S+=A-L;

                }else{

                fact1= n1_kq[l]*n1_kq[m];
                fact2= n1_kq[l]*n1_kq[m]-l_kqk1q1[l][m];
                fact3= l_kqk1q1[l][m];
               
                
                if(n1_kq[l]*n1_kq[m]-l_kqk1q1[l][m]<=1)
  
                {A=0.;}else{
                c1=1;            
                for(c = fact2; c <= fact1; c++)
                 {c1 = c1 * c;}
                 A=log(c1);
                 }
 
                if(l_kqk1q1[l][m]<=1){
                L=0.;
                }else{
                c1=1;            
                for(c = 1; c <= fact3; c++)
                 {c1 = c1 * c;}
                 L=log(c1);

                }
 
                S+=A-L;
                if(A>DBL_MAX){
                printf("Network size N is too large! Use Stirling mode...\n"); 
                exit(0);}
               }  
             //printf("l=%d m=%d n1_kq[l]=%d n1_kq[m]=%d l_kqk1q1[l][m]=%d fact1=%lf fact2=%lf fact3=%lf A=%lf B=%lf L=%lf S=%lf\n",l,m,n1_kq[l],n1_kq[m],l_kqk1q1[l][m],fact1,fact2,fact3,A,B,L,S);        
             }else{S+=0.;}
            }
           }
          }

          //for(i=0;i<maxkc+1;i++){
                       //free(n_kq[i]);}
          //free(n_kq);
          free(n1_kq);

          for(i=0;i<d+1;i++){
                       free(l_kqk1q1[i]);}
          free(l_kqk1q1);


          return(S/(double)T);
       }





/* Statistics. */

double stat(double datum, double *data, int lenght){

    int i;
    double mean, sum_deviation, Z;
       

    for(i=0; i<lenght; i++){
        mean+=data[i];
    }
    mean=mean/(double)lenght;
    //printf("Nr=%d mean=%g \n",lenght, mean);
    for(i=0; i<lenght; i++){
    sum_deviation+=pow((data[i]-mean),2);}
    //printf("sumdev=%g \n",sum_deviation);     
    Z=(datum-mean)/sqrt(sum_deviation/(double)(lenght-1));
    //printf("Z=%lf \n",Z);  
    return Z;

}





/* Main*/
int main(int argc, char* argv[]){

	int i,j,l,m,c,f,h,n,T,ci,cl,cf,cn,maxk1,maxq1,maxkc1,maxk2,maxq2,maxkc2,dim1,dim2,k_class_mode,stirling_mode,Nr,*kt1,*kt2;
        double Sigma11,Sigma22,Sigma12,Sigma21,*SigmaR11, *SigmaR22, *SigmaR12, *SigmaR21, Theta11, Theta22, Theta12, Theta21, ThetaS;	
         

        char filec[60];
        char *filename1, *filename2, *filename3, *filename4;
 

   printf("\n\n ++++++++ MEMSA v 1.0 ++++++ \n MEMSA (MEsoscopic Multiplex Structure Analysis) is a program developed to extract information from the mesoscopic structures of multiplex\n networks induced by node features.\n This program can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software \n Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed by the author in the \n hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR \n PURPOSE. If you use this program please cite \n [1] J. Iacovacci, Z. Wu and G. Bianconi, Mesoscopic structures reveal the network between the layers of multiplex datasets,\n Physical Review E 92 (4), 2015  \n\n (c) Jacopo Iacovacci (mriacovacci@hotmail.it) \n ++++++++++++++++++++++++\n\n");

printf(" Usage:\n\n ./theta_S [N] [edge_list_layer1.txt] [feature_layer1.txt] [edge_list_layer2.txt] [feature_layer2.txt] [k-classes] [stirling] [Nr] \n\n N -- number of nodes \n k-classes -- use of classes of nodes respect to the degree (0 = off, 1 = on)\n stirling -- Stirling approximation to compute the entropy (0 = off, 1 = on)\n Nr -- numbers of reshuffling to compute the similarity score\n\n");

        if(argc!=9){
		printf("Input files or values not correct. Please check\n");
		return 1;
	}


        FILE *fp;

	T=atoi(argv[1]);
        filename1 = argv[2];
        filename2 = argv[3];
        filename3 = argv[4];
        filename4 = argv[5];
        k_class_mode=atoi(argv[6]);  
        stirling_mode=atoi(argv[7]);
        Nr=atoi(argv[8]);
   
        k1=(int*)calloc(T,sizeof(int));
        kt1=(int*)calloc(T,sizeof(int));
       
        k2=(int*)calloc(T,sizeof(int));
        kt2=(int*)calloc(T,sizeof(int));


        maxk1 = 0;
        maxq1 = 0;
        maxkc1 = 0;

        maxk2 = 0;
        maxq2 = 0;
        maxkc2 = 0;

  
	if(filename1==NULL) {
		printf("\aError.\nCannot open file %s\n",filename1);}
 
        fp = fopen(filename1,"r");
	while(fscanf(fp,"%d %d\n",&i ,&j)!=EOF){
		 
		k1[i]++;
		k1[j]++;
                if(k1[i]>maxk1)maxk1=k1[i];
                if(k1[j]>maxk1)maxk1=k1[j];		 

		 
	}
	fclose(fp);
 
        knn1=(int**)calloc(T,sizeof(int*));

	for(i=0;i<T;i++){
		knn1[i]=(int*)calloc(k1[i],sizeof(int));}

        fp = fopen(filename1,"r");
	while(fscanf(fp,"%d %d\n",&i ,&j)!=EOF){
		 
		kt1[i]++;
		kt1[j]++;

                knn1[i][kt1[i]-1]=j;
	        knn1[j][kt1[j]-1]=i;		 
	}
	fclose(fp);

        free(kt1);



	if(filename3==NULL) {
		printf("\aError.\nCannot open file %s\n",filename3);}

        fp = fopen(filename3,"r");
	while(fscanf(fp,"%d %d\n",&i ,&j)!=EOF){
		 
		k2[i]++;
		k2[j]++;
                if(k2[i]>maxk2)maxk2=k2[i];
                if(k2[j]>maxk2)maxk2=k2[j];		 
	 
	}
	fclose(fp);


        knn2=(int**)calloc(T,sizeof(int*));

	for(i=0;i<T;i++){
		knn2[i]=(int*)calloc(k2[i],sizeof(int));}

        fp = fopen(filename3,"r");
	while(fscanf(fp,"%d %d\n",&i ,&j)!=EOF){
		 
		kt2[i]++;
		kt2[j]++;

                knn2[i][kt2[i]-1]=j;
	        knn2[j][kt2[j]-1]=i;		 
	}
	fclose(fp);

        free(kt2);


	q1=(int*)calloc(T,sizeof(int));
        q2=(int*)calloc(T,sizeof(int));


        if(filename2==NULL) {
		printf("\aError.\nCannot open file %s\n",filename2);}

        fp=fopen(filename2, "r");
	while(fscanf(fp, "%d %d",&i,&cl)!=EOF){
		q1[i]=cl;
	        if(cl>maxq1)maxq1=cl;}
	fclose(fp);

        if(filename4==NULL) {
		printf("\aError.\nCannot open file %s\n",filename4);}
  
        fp=fopen(filename4, "r");
	while(fscanf(fp, "%d %d",&i,&cl)!=EOF){
		q2[i]=cl;
	        if(cl>maxq2)maxq2=cl;}
	fclose(fp);


        kc1=(int*)calloc(T,sizeof(int));
        kc2=(int*)calloc(T,sizeof(int));
 
        if(k_class_mode==1){
        
        maxkc1=classes(k1, maxk1, T, kc1);
        maxkc2=classes(k2, maxk2, T, kc2);

        }else{

        for(i=0;i<T;i++){
        
        kc1[i]=k1[i];
        kc2[i]=k2[i];}
        
        maxkc1=maxk1;
        maxkc2=maxk2;

        }

        dim1 = maxkc1*maxq1;
        dim2 = maxkc2*maxq2;
        printf("\n", dim1, dim2);
        printf("dimensionality of layer 1 = %d  \ndimensionality of layer 2 = %d\n", dim1, dim2);
        printf("\n", dim1, dim2);
        //printf("#classes_layer1 x #_communities_layer1=%d  #classes_layer2 x #communities_layer2=%d\n", dim1, dim2);
        //printf("Ecco maxk1=%d maxk2=%d maxkc1=%d maxkc2=%d\n",maxk1, maxk2,maxkc1, maxkc2);

        
        Sigma11=entropy(T, knn1, k1, kc1, q1, maxk1, maxkc1, maxq1, stirling_mode);
        Sigma22=entropy(T, knn2, k2, kc2, q2, maxk2, maxkc2, maxq2, stirling_mode); 
        Sigma12=entropy(T, knn1, k1, kc1, q2, maxk1, maxkc1, maxq2, stirling_mode);
        Sigma21=entropy(T, knn2, k2, kc2, q1, maxk2, maxkc2, maxq1, stirling_mode); 
  
        //printf("%g %g % g %g \n",Sigma11, Sigma22, Sigma12, Sigma21); 
        //return 0;

        SigmaR11=(double*)calloc(Nr,sizeof(double));
        SigmaR22=(double*)calloc(Nr,sizeof(double));
        SigmaR12=(double*)calloc(Nr,sizeof(double));
        SigmaR21=(double*)calloc(Nr,sizeof(double));



        for(n=0;n<Nr;n++){

        shuffle(q1, T);
        shuffle(q2, T);

        SigmaR11[n]=entropy(T, knn1, k1, kc1, q1, maxk1, maxkc1, maxq1, stirling_mode);
        SigmaR22[n]=entropy(T, knn2, k2, kc2, q2, maxk2, maxkc2, maxq2, stirling_mode);   
        SigmaR12[n]=entropy(T, knn1, k1, kc1, q2, maxk1, maxkc1, maxq2, stirling_mode);
        SigmaR21[n]=entropy(T, knn2, k2, kc2, q1, maxk2, maxkc2, maxq1, stirling_mode);   

        //printf("%g %g % g %g \n",SigmaR11[n], SigmaR22[n], SigmaR12[n], SigmaR21[n]);
        }
  
        Theta11=stat(Sigma11, SigmaR11, Nr);
        Theta22=stat(Sigma22, SigmaR22, Nr);

        Theta12=stat(Sigma12, SigmaR12, Nr);
        Theta21=stat(Sigma21, SigmaR21, Nr);
 

        ThetaS=(Theta12/Theta11+Theta21/Theta22)/2;
       
        //printf("Theta11 -- Theta22 -- Theta12 -- Theta21 -- ThetaS \n");
        //printf("%g   %g   %g   %g   %g\n",Theta11, Theta22, Theta12, Theta21, ThetaS);

        
        printf("Theta(1,1) = %g \n",Theta11);
        printf("Theta(1,2) = %g \n",Theta12);
        printf("Theta(2,1) = %g \n",Theta21);
        printf("Theta(2,2) = %g \n\n",Theta22);
        printf("ThetaS score = %g \n\n",ThetaS);
        //printf("%g   %g   %g   %g   %g\n",Theta11, Theta22, Theta12, Theta21, ThetaS);
    
        for(i=0;i<T;i++){
                       free(knn1[i]);
                       free(knn2[i]);}
        free(knn1);
        free(knn2);
  
        free(k1);
        free(k2);
        free(q1);
        free(q2);
        free(kc1);
        free(kc2);
        free(SigmaR11);
        free(SigmaR22);
        free(SigmaR12);
        free(SigmaR21);

      return 0;
  }




