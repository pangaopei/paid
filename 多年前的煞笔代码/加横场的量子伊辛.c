#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define N 16 //N^2个粒子
int hm = 10; //h最大步数
double dh = 0.04; //h步长
int ustep = 80000; //模拟次数
double eps = 1e-8;
int s[N][N]; //自旋
int m0; //总自旋
double h=0.8;
double P; //哈密顿量
double J = 44;
double d_tau = 0.01;
double H;
// beta = N * d_tau

//初始化 
void init(double h){
    m0 = 0;
    int i,j; 
  
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++){
			s[i][j] = (rand() % 2) ? 1 : -1;
    		m0 += s[i][j];
		}
	H = 0;
	for (i = 0; i < N; i++)
		for( j=0;j<N;j++)
    		H -= J*s[i][j] * s[(i + 1) % N][j];
}

//改变(i,j)处自旋方向的改变
double deltaP(int i, int j,double h){
    return -2 * s[i][j] *((s[(i + 1) % N][j] +s[(i - 1 + N) % N][j] )*d_tau * J-0.5*log(tanh(h*d_tau))*(s[i][(j + 1) % N] + s[i][(j - 1 + N) % N]));
}

void Metropolis(double h){
    double dP = 0; //改变量
    double pr; // [0,1)随机数
    int i,p,q; 

    for (i = 0; i < 5; i++){
    	for(p=0;p<N;p++){
    		for(q=0;q<N;q++){
				dP = deltaP(p, q,h);
				pr = (double)(rand()) / (double)(RAND_MAX);
				if (dP > eps || pr  - exp(dP) < - eps){
					s[p][q] = - s[p][q];
					m0 += 2 * s[p][q];
					H -=2*J * s[p][q] * (s[(p + 1) % N][q] + s[(p-1+N)%N][q]);
				}	
			}
		} 
    }
}


int main(){
	double m_bar[hm];
	double m2[hm];
	double d[hm];
	double H_bar[hm];
	double h2[hm];
	double dd[hm];
    FILE *fp;
		int i;
 		char name[30];
 		printf("请输入文件名  :");
 		scanf("%s",name);
		fp = fopen(name, "w+");
int hstep=0;
		for (hstep = 0; hstep < hm; hstep++, h += dh){
        m_bar[hstep] = 0;
		m2[hstep]=0;
		d[hstep]=0;
		dd[hstep]=0;
		h2[hstep]=0;
		H_bar[hstep] = 0;
		    srand((int)time(0));
        init(h*J);
            //平衡 
            int step;
            for (step = 0; step < 80000; step++) Metropolis(h*J);

            //模拟ustep次
            for (step = 0; step < ustep; step++){
                Metropolis(h*J);
                m_bar[hstep] +=fabs((double)(m0))/(double)(N)/(double)(N);
				m2[hstep]+=(double)(m0)*(double)(m0)/(double)(N)/(double)(N)/(double)(N)/(double)(N);
				H_bar[hstep] +=(double)(H)/(double)(N)/(double)(N);
				h2[hstep] +=(double)(H)*(double)(H)/(double)(N)/(double)(N)/(double)(N)/(double)(N);
            }
        m_bar[hstep] /= (double)(ustep) ;
		m2[hstep] /= (double)(ustep);
		d[hstep]=m2[hstep]-m_bar[hstep]*m_bar[hstep];
		H_bar[hstep] /= (double)(ustep)*(double)(J);
		h2[hstep] /= (double)(ustep)*(double)(J);
		dd[hstep]=h2[hstep]-H_bar[hstep]*H_bar[hstep];
        fprintf(fp,"%lf  %lf  %lf  %lf  %lf\n",h,fabs(m_bar[hstep]),d[hstep],H_bar[hstep],dd[hstep]);
        
        printf("%.2lf%%\th=%.3lf\tm=%.4lf\td=%.4lf\nH=%.4lf\tdd=%.4lf\n",(double)(hstep+1) / (double)(hm) * 100.0, h, fabs(m_bar[hstep]),d[hstep],H_bar[hstep],dd[hstep]);
    }
    fclose(fp);
}

