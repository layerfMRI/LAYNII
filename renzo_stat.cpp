

#include "renzo_stat.h"  
//#include <math.h>

  
  double ren_average(double arr[], int size) {
  int i;
  double sum = 0;       
  double avg;          

   for (i = 0; i < size; ++i) {
      sum += arr[i];
   }
   avg = double(sum) / size;

   return avg;
}


double ren_stdev(double arr[], int size) {
   
  int i;
  double sum = 0;    
  double mean ;    
 
  
  mean = ren_average(arr, size) ; 
  
  for (i = 0; i < size; ++i) {
      sum += (arr[i]-mean)*(arr[i]-mean)/((double)size-1);
  }
	
  return sqrt( sum ) ; 
  	
}


double ren_correl(double arr1[], double arr2[], int size) {
	
	int i;
	double sum1 = 0 ;
	double sum2 = 0 ; 
	double sum3 = 0 ; 	
	double mean1 = ren_average(arr1, size) ; 
	double mean2 = ren_average(arr2, size) ; 	
	
	  for (i = 0; i < size; ++i) {
       sum1 += (arr1[i]-mean1)*(arr2[i]-mean2);
       sum2 += (arr1[i]-mean1)*(arr1[i]-mean1);
       sum3 += (arr2[i]-mean2)*(arr2[i]-mean2);
      }
	
	
	return sum1/sqrt(sum2*sum3); 
	
}


double ren_skew(double arr[], int size) {
  int i;
  double sum1 = 0;   
  double sum2 = 0;   
  double mean = ren_average(arr, size) ; 
  
   for (i = 0; i < size; ++i) {
       sum1 += (arr[i]-mean)*(arr[i]-mean)*(arr[i]-mean);
       sum2 += (arr[i]-mean)*(arr[i]-mean);
      }
	
  return (1/((double)size) * sum1 )/( pow ( 1/((double)size-1) * sum2 , 1.5 )  ) ;  
  
}

double ren_kurt(double arr[], int size)  {
  int i;
  double sum1 = 0;    
  double sum2 = 0; 
  double mean = ren_average(arr, size) ; 
    
     for (i = 0; i < size; ++i) {
       sum1 += (arr[i]-mean)*(arr[i]-mean)*(arr[i]-mean)*(arr[i]-mean)/((double)size);
       sum2 += (arr[i]-mean)*(arr[i]-mean)/((double)size);
      }
  
  return sum1/(sum2*sum2)-3 ; 
}

double ren_autocor(double arr[], int size)  {
  int i;
  double sum1 = 0;    
  double sum2 = 0; 
  double mean = ren_average(arr, size) ; 
  
      for (i = 1; i < size; ++i) {
        sum1 += (arr[i]-mean)*(arr[i-1]-mean);
      }
  
      for (i = 0; i < size; ++i) {
        sum2 += (arr[i]-mean)*(arr[i]-mean);
      }
  
  return sum1/sum2 ; 
}
