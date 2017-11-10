#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#define PI 3.14159265359
 /*-------------------
 *
 *
 *      Libraries
 *
 *
 * -------------------*/
#define alpha_m
#define ma 0.01
#define mc 0.01 
#define R0 (1-(double)(alpha_m))log2(1+(double)(gam_m)) //the throughput in the hop
#define gam_m //gamma_m
#define n0  //noise
#define delta 1 //variance of gaussian distribution
#define epsilon 0.5 //parameter of estimated channel
double bessel1k(double n);
double uniform(void);
double Box(double s, double r1, double r2, double *x,double *y);
  int i;
  double s;




/*----------------------------------------------------*
 *                                                    *
 *                                                    *
 *      Rayleigh fading random variables              *
 *                                                    *
 *                                                    *
 *----------------------------------------------------*/
double rayleigh(void)
{
  double x,y,h;
  double r1,r2;
  r1=(double)random()+1.0/((double) RAND_MAX+1.0);
  r2=(double)random()+1.0/((double) RAND_MAX+1.0);
  x = s * sqrt( -2.0 * log (r1)) * cos (2.0 * M_PI * r2 );
  y = s * sqrt( -2.0 * log (r1)) * cos (2.0 * M_PI * r2 );
  h = pow(x,2.0)+pow(y,2.0);
  return h;
}



/*----------------------------------------------------*
 *                                                    *
 *                                                    *
 *     log-normal fading random variables             *
 *                                                    *
 *                                                    *
 *----------------------------------------------------*/

int main(){
  srandom((unsigned)time(NULL));
  /*------------------------------------------
   *
   *  Define parameters
   *
   *------------------------------------------*/
  double d1,d2,m,delta_r,delta_d,gamma_o,eta,SNR;
  double a,b,c,d,u;
  double tao[200],tao1[200],tao2[200],tao3[200],tao4[200],tao5[200],tao6[200],tao7[200],tao8[200],tao9[200];
  double count[200],count1[200];
  double p_Ps,p_gamma_o,p_delta_r,p_delta_d,p_Pr;
  double p_gamma_di,alpha;
  double r_h,p_h,h_i_2,r_g,p_g,g_i_2,h_m,g_m;
  double tao_i_s,tao_i_d,tao_i_sd,E_i_s,E_i_d,E_i_sd,alpha_i_s,alpha_i_d,alpha_i_sd,I_o,E_i_new,alpha_i_new,tao_i_new;
  double T,en;
  double alpha_P;
  double c_counter;
  int N,i,counter,Pr,Ps;
  double pr_o_s,pr_o_d,pr_o_sd;
  alpha=0.5;
  counter=0;
  //Pr = 0;//db
  Ps = 30;
  gamma_o = 10;
  d1 = 200;
  d2 = 50;
  alpha_P = 0.05;
  m = 3.0;
  delta_r = -100.0;
  delta_d = -100.0;
  eta = 0.5;
  en = 0;
  /*------------------------------------------
   *
   *
   * Throughput efficiency
   *
   *
   *------------------------------------------*/

  //analytical
  
  p_gamma_o = pow(10,0.1*gamma_o);
  p_delta_r = pow(10,-3)*pow(10,0.1*delta_r);
  //p_delta_r = delta_r;
  //p_delta_d = delta_d;
  p_delta_d = pow(10,-3)*pow(10,0.1*delta_d);
  
 for(Pr = -45; Pr <= -25; Pr++)
  {
    //d2 = 100 - d1;
    //p_Ps = pow(10,0.1*(double)(Ps));
    p_Ps = pow(10,-3)*pow(10,0.1*(double)(Ps));//Ps has been changed.
    p_Pr = pow(10,-3)*pow(10,0.1*(double)(Pr));
    a = p_Ps * pow(d2,m) * p_delta_d * p_gamma_o;
    b = pow(d1,m) * pow(d2,m) * p_delta_r * p_delta_d * p_gamma_o;
    c = p_Ps * p_Pr;
    d = p_Pr * pow(d1,m) * p_delta_r * p_gamma_o;
    u = sqrt(4 * ( a * d + b * c)/ pow(c,2));
    //SD-WPT
    tao6[(int)(Pr+45)] = exp(-(a+d)/c) * u * bessel1k(u)/2 / (1 + pow(p_Pr,2) *pow(d1,m)*pow(d2,m)/(4*(pow(eta,2))*pow(p_Ps,2)));
    //tao8[(int)(Pr+45)] = exp(-(a+d)/c) / 2 / (1 + p_Pr *pow(d1,m)*pow(d2,m)/(2*(eta*pow(d2,m)+pow(d1,m)*eta)*p_Ps));
    //S-WPT
    tao[(int)(Pr+45)] = exp(-(a+d)/c) * u * bessel1k(u)/ 2 / (1 + p_Pr *pow(d1,m)/(2*eta*p_Ps));
    tao2[(int)(Pr+45)] = exp(-(a+d)/c) /2/( 1+ p_Pr * pow(d1,m)/(2 * eta * p_Ps));
    //D-WPT
    tao3[(int)(Pr+45)] = exp(-(a+d)/c) * u * bessel1k(u)/ 2 / (1 + p_Pr *pow(d2,m)/(2*eta*p_Ps));
    tao5[(int)(Pr+45)] = exp(-(a+d)/c) / 2 / (1 + p_Pr *pow(d2,m)/(2*eta*p_Ps));
    //S-WPT optimal
    pr_o_s = 10*log10(pow(10,3)*pow(d1,-m)*(0.5*(pow(d1,m)*pow(d2,m)*p_gamma_o*p_delta_d+sqrt(p_gamma_o*pow(d1,m)*pow(d2,m)*p_delta_d*(8*eta*p_Ps+pow(d1,m)*pow(d2,m)*p_gamma_o*p_delta_d)))));
    pr_o_d = 10*log10(pow(10,3)*0.5*(pow(d2,m) * p_gamma_o * p_delta_d+sqrt(p_gamma_o * p_delta_d *(8*eta*p_Ps+pow(d2,2*m)*p_gamma_o*p_delta_d))));
    pr_o_sd = 10*log10(pow(10,3)*0.5*pow(d1,-m)*(pow(d1,m)*pow(d2,m)*p_gamma_o*p_delta_d+sqrt(pow(d1,m)*p_gamma_o*p_delta_d*(8*(1-alpha)*pow(d1,m)*eta*p_Ps+8*alpha*pow(d2,m)*eta*p_Ps+pow(d1,m)*pow(d2,2*m)*p_gamma_o*p_delta_d))));
      
    
     printf("Pr    u    uK(u)    tao_s    tao_sa   tao_d   tao_da   tao_sd   tao_sda\n");
     printf("%d %f %f %f %f %f %f %f %f\n",Pr,u,u*bessel1k(u),tao[(int)(Pr+45)],tao2[(int)(Pr+45)],tao3[(int)(Pr+45)],tao5[(int)(Pr+45)],tao6[(int)(Pr+45)],tao8[(int)(Pr+45)]);
     //Optimal Pr

     printf("S-WPT    D-WPT    SD-WPT\n");
     printf("%f %f %f\n",pr_o_s,pr_o_d,pr_o_sd);
  //printf("P_PS  P_gamma_o p_delta_d p_delta_r u*y1(u) p_pr\n");
  //printf("%f %f %f %f %f %f\n",p_Ps,p_gamma_o,log(p_delta_d),log(p_delta_r),u*y1(u)*(exp(-(a+d)/c)),0.1*(double)(Pr));
  //printf("%f\n",-(a+d)/c);
  }
  
  //simulation
  double h_i[100],g_i[100];
  double ene[200];
  int n;
  c_counter=0;
  N = 100000;

  T = 10000;//simulation time
  
  
  for(Pr = -45;Pr <= -25; Pr++){
    tao_i_s= 0.0;
    tao_i_d = 0.0;
    tao_i_sd = 0.0;
    tao_i_new = 0.0;
    counter=0;
    SNR=0;
    //d2 = 100-d1;
    en=0;
    E_i_s = 0;
    E_i_d = 0;
    E_i_sd = 0;
    E_i_new = 0;
    h_m=0;//initial parameter of memory
    g_m=0;//
  for (i=0;i<N;i++)
    {      
  double uniform(void);
   p_Pr = pow(10,-3) * pow(10,0.1*(double)(Pr));
  // p_Pr = pow(10,0.1 * 5) / sqrt(0.5 * PI);// 5db/exceptional number of rayleigh fading
  p_Ps = pow(10,-3) * pow(10,0.1*(double)(Ps));
  // p_Ps = pow(10,0.1*(double)(Ps));
  Box(1,uniform(),uniform(),&r_h,&p_h);
  Box(1,uniform(),uniform(),&r_g,&p_g);
  h_i_2 = (pow(r_h,2) +pow(p_h,2))/2;
  //h_i_2 = (1-pow(epsilon,2))*(pow(r_h,2) + pow(p_h,2))/2 + pow(epsilon,2) * h_m +
  //   2 * epsilon * sqrt(1-pow(epsilon,2))* sqrt(h_m)*sqrt((pow(r_h,2)+pow(p_h,2))/2);//rayleigh fading parameter square
  g_i_2 =(pow(r_g,2) + pow(p_g,2))/2;//rayleigh fading parameter square
  //g_i_2 = (1-pow(epsilon,2))*(pow(r_g,2) + pow(p_g,2))/2 + pow(epsilon,2) * g_m +
  //  2 * epsilon * sqrt(1-pow(epsilon,2))* sqrt(g_m)*sqrt((pow(r_g,2)+pow(p_g,2))/2);
  h_m = h_i_2;//memory latest statement of fading parameter
  g_m = g_i_2;
  // printf("%f %f\n",h_i_2,g_i_2);



  //log-normal miu=0 sigma=1
  // h_i_2 = pow(exp(r_h),2);
  //g_i_2 = pow(exp(p_h),2);


  p_gamma_di = p_Ps * p_Pr * h_i_2 * g_i_2 /(p_Pr * g_i_2 * pow(d1,m) * p_delta_r + p_delta_d * pow(d2,m) *(p_Ps *h_i_2 + pow(d1,m)* p_delta_r));
  SNR += p_gamma_di;

   if(p_gamma_di < p_gamma_o)
    {
    I_o = 1;
    //c_counter += 1;
    }
  else
    {
    I_o = 0;

    }
   
  //after T time
   if(E_i_s >= p_Pr/2  )
    {
       alpha_i_s = 0;
       // E_i_s = E_i_s - p_Pr/2 ;
       E_i_s = 0;
    //counter=0;
     }
     else
    {
      //S-WPT
      E_i_s = E_i_s + eta * p_Ps * (h_i_2)/pow(d1,m); 
      
      alpha_i_s = 1;
      counter++;
      }
  
   if (alpha_i_s == 0 && I_o == 0){

     tao_i_s += 0.5;
     //counter++;
   }

   
   //D-WPT
   if(E_i_d >= p_Pr/2  )
     {
       alpha_i_d = 0;
       E_i_d = E_i_d - p_Pr/2 ;
       //counter=0;
     }
     else
    {
      //D-WPT
       E_i_d = E_i_d +eta * p_Ps * (g_i_2)/pow(d2,m);
      
      alpha_i_d = 1;
      counter++;
      }
  
   if (alpha_i_d == 0 && I_o == 0){

     tao_i_d += 0.5;
     //counter++;
   }

   //SD-WPT
   if(E_i_sd >= p_Pr/2  )
    {
       alpha_i_sd = 0;
       // E_i_sd = E_i_sd - p_Pr/2 ;
       E_i_sd = 0;
    //counter=0;
     }
     else
    {
      //SD-WPT
       E_i_sd = E_i_sd + eta * p_Ps * (h_i_2)/pow(d1,m) + eta * p_Ps * (g_i_2)/pow(d2,m);

      alpha_i_sd = 1;
      counter++;
      }
  
   if (alpha_i_sd == 0 && I_o == 0){

     tao_i_sd += 0.5;
     //counter++;
   }



   //new-WPT
   if(E_i_new >= p_Pr/2 )
     {
       alpha_i_new = 0;
       E_i_new = E_i_new - p_Pr/2;
       //counter=0;
     }
     else
    {
      //new-WPT
      if((h_i_2/pow(d1,m)) > (g_i_2/pow(d2,m)))
	{
       E_i_new = E_i_new +eta * p_Ps * (h_i_2)/pow(d1,m);
       alpha_i_new = 1;
      /* counter++; */
	}
      else{

	E_i_new = E_i_new + eta * p_Ps * (g_i_2)/pow(d2,m);
	alpha_i_new = 1;
      }
    }

   if(alpha_i_new ==0 && I_o == 0){

     tao_i_new += 0.5;
   }
  

   //printf("%d %d\n",i,counter);

    }
  SNR = SNR/N;
  printf("Ps:%d SNR:%f Pr:%d\n",Ps,10*log(SNR),Pr);
  c_counter = (double)(counter)/(double)(N);
  tao_i_s =  tao_i_s / (double)(N);
  tao_i_d = tao_i_d / (double)(N);
  tao_i_sd = tao_i_sd / (double)(N);
  tao_i_new = tao_i_new / (double)(N);
  tao1[(int)(Pr+45)] = tao_i_s;
  tao4[(int)(Pr+45)] = tao_i_d;
  tao7[(int)(Pr+45)] = tao_i_sd;
  tao9[(int)(Pr+45)] = tao_i_new;
  count[(int)(Pr+45)] = c_counter;
  // ene[(int)(d1)]= en/counter;
  printf("Pr      S-WPT     D-WPT     SD-WPT      new\n");
  printf("%d %f %f %f %f\n",(int)(Pr),tao1[(int)(Pr+45)],tao4[(int)(Pr+45)],tao7[(int)(Pr+45)],tao9[(int)(Pr+45)]);
  }

  

  //result
  FILE *file;
  file=fopen("throughputefficiency_s.txt","w");
  for(Pr = -45; Pr <= -25; Pr++)
    {
      fprintf(file,"%d %f\n",(int)(Pr),tao[(int)(Pr+45)]);
    }
  fclose(file);
    
  file=fopen("throughputefficiency1_s.txt","w");
  for(Pr = -45; Pr <= -25; Pr++)
    {
      fprintf(file,"%d %f\n",(int)(Pr),tao1[(int)(Pr+45)]);
    }
  fclose(file);
  
  file=fopen("throughputefficiency2_s.txt","w");
  for(Pr = -45;Pr <= -25;Pr++)
    {
      fprintf(file,"%d %f\n",(int)(Pr),tao2[(int)(Pr+45)]);
    }
  fclose(file);
    file=fopen("throughputefficiency_d.txt","w");
  for(Pr = -45; Pr <= -25; Pr++)
    {
      fprintf(file,"%d %f\n",(int)(Pr),tao3[(int)(Pr+45)]);
    }
  fclose(file);
    
  file=fopen("throughputefficiency1_d.txt","w");
  for(Pr = -45; Pr <= -25; Pr++)
    {
      fprintf(file,"%d %f\n",(int)(Pr),tao4[(int)(Pr+45)]);
    }
  fclose(file);
  
  file=fopen("throughputefficiency2_d.txt","w");
  for(Pr = -45;Pr <= -25;Pr++)
    {
      fprintf(file,"%d %f\n",(int)(Pr),tao5[(int)(Pr+45)]);
    }
  fclose(file);
    file=fopen("throughputefficiency_sd.txt","w");
  for(Pr = -45; Pr <= -25; Pr++)
    {
      fprintf(file,"%d %f\n",(int)(Pr),tao6[(int)(Pr+45)]);
    }
  fclose(file);
    
  file=fopen("throughputefficiency1_sd.txt","w");
  for(Pr = -45; Pr <= -25; Pr++)
    {
      fprintf(file,"%d %f\n",(int)(Pr),tao7[(int)(Pr+45)]);
    }
  fclose(file);
  
  file=fopen("throughputefficiency2_sd.txt","w");
  for(Pr = -45;Pr <= -25;Pr++)
    {
      fprintf(file,"%d %f\n",(int)(Pr),tao8[(int)(Pr+45)]);
    }
  fclose(file);
  file=fopen("throughputefficiency_new.txt","w");
  for(Pr = -45;Pr <= -25;Pr++)
    {
      fprintf(file,"%d %f\n",(int)(Pr),tao9[(int)(Pr+45)]);
    }
  fclose(file);
  /* file=fopen("energy.txt","w");
  for(d1 = 2; d1 <= 48; d1++)
    {
      fprintf(file,"%d %f\n",(int)(d1),ene[(int)(d1)]);
    }
    fclose(file);*/
  /*
  file=fopen("randomh.txt","w");
  for(n = 0; n < 100; n++)
    {
      fprintf(file,"%f %f\n",(double)(n)/10,(double)(h_i[n])/(double)(N));
    }
  fclose(file);
  file=fopen("randomg.txt","w");
  for(n = 0; n<100;n++)
    {
      fprintf(file,"%f %f\n",(double)(n)/10,(double)(g_i[n])/(double)(N));
    }
    fclose(file);*/

  return 0;
}



/*Uniform pseudorandom number used random function */
/*(0,1]*/
double uniform(void){
  return((double)random() + 1.0) / ((double)RAND_MAX + 1.0);
}

/*Box-muller's method'*/
/*r1,r2は一様乱数,x,yは正規乱数の出力先*/
/*平均値は0,分散がs^2*/
double Box(double s,double r1,double r2,double *x,double *y){
  *x = s * sqrt(-2.0*log(r1)) * cos(2.0*M_PI*r2);
  *y = s * sqrt(-2.0*log(r1)) * sin(2.0*M_PI*r2);
  return 0.0;
}


/* Modified Bessel function of the second kind*/
double besseli1(double n)
{
  double ax,ans;
  double y;

  if((ax = fabs(n)) < 3.75){
    y = n/3.75;
    y*=y;
    ans = ax * (0.5 + y * (0.87890594 + y * (0.51498868 + y * (0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
  }
  else{
    y = 3.75/ax;
    ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
    ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y *ans))));
    ans *= (exp(ax)/sqrt(ax));
  }
  return n < 0.0 ? -ans : ans;
}

double bessel1k(double n){
  double besseli1(double n);
  double y,ans;

  if(n <= 2.0){
    y = n * n / 4.0;
    ans = (log ( n / 2.0 ) * besseli1(n)) + (1.0 / n) * (1.0 + y * (0.15443144 + y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1 + y * (-0.110404e-2 + y * (-0.4686e-4)))))));
  }
  else {
      y = 2.0 / n;
      ans = (exp(-n) / sqrt(n)) * (1.25331414 + y * (0.23498619 + y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2 + y * (0.325614e-2 + y * (-0.68245e-3)))))));
  }
  return ans;
}
