#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N 32 //N^2个粒子
////int Tm = 80; //T最大步数
double dT = 0.05; //T步长
int ustep = 200000; //达到热平衡后的模拟次数
double eps = 1e-8;
int s[N][N]; //自旋
int m0; //总自旋
double H; //哈密顿量
double h = 0.0; //外场项
double T = 2; // (J / kB)
//J = 1 
int k;

//初始化 
void init(){
    m0 = 0;
    int i,j; 
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++){
                s[i][j] = (rand() % 2) ? 1 : -1;
                m0 += s[i][j];
            }
    H = 0;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
                H -= s[i][j] * (s[(i + 1) % N][j] + s[i][(j + 1) % N] + h);
}

//改变(i,j)处自旋方向能量的改变
double deltaH(int i, int j){
    return 2 * s[i][j] * (s[(i + 1) % N][j] + s[i][(j + 1) % N] +
                             s[(i - 1 + N) % N][j] + s[i][(j - 1 + N) % N] + h);
}
void new(int i,int j,int k){
	double pr; // [0,1)随机数
	if(s[(i+1)%N][j]==k){
		pr = (double)(rand()) / (double)(RAND_MAX);
		if (pr  -1.0+ exp(-2/T) < - eps){
            	s[(i+1)%N][j]*=-1;
				new((i+1)%N,j,k);	
        }
	}
	if(s[(i-1+N)%N][j]==k){
		pr = (double)(rand()) / (double)(RAND_MAX);
		if (pr  -1.0+ exp(-2/T) < - eps){
            	s[(i-1+N)%N][j]*=-1;
				new((i-1+N)%N,j,k);	
        }
	}
	if(s[i][(j+1)%N]==k){
		pr = (double)(rand()) / (double)(RAND_MAX);
		if (pr  -1.0+ exp(-2/T) < - eps){
            	s[i][(j+1)%N]*=-1;
				new(i,(j+1)%N,k);	
        }
	}
	if(s[i][(j-1+N)%N]==k){
		pr = (double)(rand()) / (double)(RAND_MAX);
		if (pr  -1.0+ exp(-2/T) < - eps){
            	s[i][(j-1+N)%N]*=-1;
				new(i,(j-1+N)%N,k);	
        }
	}
}

void Metro(){
 //   for (i = 0; i < N*N; i++){
        int p = rand() % N, q = rand() % N;
		k=s[p][q];
		s[p][q]=-k;
		H+=deltaH(p,q);
		m0 += 2 * s[p][q];
		new(p,q,k);
  //  }
}

int main(){
  //  double H_bar[Tm], H_square_bar[Tm]; //能量均值和能量平方的均值
  //  double C[Tm]; // 热容*玻尔兹曼常数
    double m_bar[ustep]; //磁矩和磁矩平方的均值
  //  double Chi[Tm]; //磁化率

    FILE *fp;
		int i;
		double t[250];
		double tau=0,kk=0,fang=0;
 		char name[30];
 		printf("请输入文件名  :");
 		scanf("%s",name);
		fp = fopen(name, "w+");
	//int Tstep;
   // for (Tstep = 1; Tstep <= Tm; Tstep++, T += dT){
      //  H_bar[Tstep] = 0;
      //  H_square_bar[Tstep] = 0;
   //     m_bar[ustep] = 0;
   //     m_square_bar[ustep] = 0;
	srand((int)time(0));
        init();
            //热平衡
            int step;
            for (step = 0; step < 50000; step++) Metro();

            //模拟ustep次
            for (step = 0; step < ustep; step++){
                Metro();
                m0 = 0;
   				int i,j; 
				for (i = 0; i < N; i++)
 				for (j = 0; j < N; j++){
					m0 += s[i][j];
            	}
   				H = 0;
    			for (i = 0; i < N; i++)
        		for (j = 0; j < N; j++)
                	H -= s[i][j] * (s[(i + 1) % N][j] + s[i][(j + 1) % N] + h);
              //  H_bar[Tstep] += H;
              //  H_square_bar[Tstep] += H * H;
                m_bar[step] =fabs((double)(m0));
                kk += fabs((double)(m0));
                fang += m0*m0;
            }
            kk /= (double)(ustep) * (double)(N * N);
            fang /=(double)(ustep) * (double)(N * N)* (double)(N * N);
            for(i=1;i<251;i++){
            	tau=0;
            	for (step = 0; step < ustep; step++){
            		tau+=m_bar[step]*m_bar[(step+i)%ustep];
            	}
            	tau /= (double)(ustep) * (double)(N * N)* (double)(N * N);
							t[i-1]=(tau-kk*kk)/(fang-kk*kk);
            
			       // H_bar[Tstep] /= (double)(ustep);
   						   //  H_square_bar[Tstep] /= (double)(ustep);
 							  //     m_bar[Tstep] /= (double)(ustep) * (double)(N * N);
   					  //   m_square_bar[Tstep] /= (double)(ustep) * (double)(N * N) * (double)(N * N);
   						    // C[Tstep] = (H_square_bar[Tstep] - H_bar[Tstep] * H_bar[Tstep]) / T / T;
     					  // Chi[Tstep] = (m_square_bar[Tstep] - m_bar[Tstep] * m_bar[Tstep]) / T * 1000.0;

       				 fprintf(fp,"%d  %lf\n",i,t[i-1]);
        
     				   printf("i=%d\t Auto=%.4lf\n",i,t[i-1]);
             }
   // }
    fclose(fp);
}

