#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N 32 //N^2个粒子
int Tm = 60; //T最大步数
double dT = 0.04; //T步长
int ustep = 500000; //达到热平衡后的模拟次数
double eps = 1e-8;
int s[N][N]; //自旋
int m0; //总自旋
double H; //哈密顿量
double h = 0.0; //外场项
double T = 0.5; //T (J / kB)
//J = 1 

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

void Metropolis(){
    double dH = 0; //能量的改变量
    double pr; // [0,1)随机数
	int i;
    for (i = 0; i < N*N; i++){
        int p = rand() % N, q = rand() % N;
        dH = deltaH(p, q);
        pr = (double)(rand()) / (double)(RAND_MAX);
        if (dH < eps || pr  - exp(- dH / T) < - eps){
            s[p][q] = - s[p][q];
            H += dH;
            m0 += 2 * s[p][q];
        }
    }
}


int main(){
    double H_bar[Tm], H_square_bar[Tm]; //能量均值和能量平方的均值
    double C[Tm]; // 热容*玻尔兹曼常数
    double m_bar[Tm], m_square_bar[Tm]; //磁矩和磁矩平方的均值
    double Chi[Tm]; //磁化率

    FILE *fp;
		int i;
 		char name[30];
 		printf("请输入文件名  :");
 		scanf("%s",name);
		fp = fopen(name, "w+");
	int Tstep;
    for (Tstep = 1; Tstep <= Tm; Tstep++, T += dT){
        H_bar[Tstep] = 0;
        H_square_bar[Tstep] = 0;
        m_bar[Tstep] = 0;
        m_square_bar[Tstep] = 0;

        init();
            //热平衡
            int step;
            for (step = 0; step < 500000; step++) Metropolis();

            //模拟ustep次
            for (step = 0; step < ustep; step++){
                Metropolis();
                H_bar[Tstep] += H;
                H_square_bar[Tstep] += H * H;
                m_bar[Tstep] += (double)(m0);
                m_square_bar[Tstep] += (double)(m0) * (double)(m0);
            }

        H_bar[Tstep] /= (double)(ustep);
        H_square_bar[Tstep] /= (double)(ustep);
        m_bar[Tstep] /= (double)(ustep) * (double)(N * N);
        m_square_bar[Tstep] /= (double)(ustep) * (double)(N * N) * (double)(N * N);
        C[Tstep] = (H_square_bar[Tstep] - H_bar[Tstep] * H_bar[Tstep]) / T / T;
        Chi[Tstep] = (m_square_bar[Tstep] - m_bar[Tstep] * m_bar[Tstep]) / T * 1000.0;

        fprintf(fp,"%lf  %lf  %lf  %lf  %lf\n",T,H_bar[Tstep],C[Tstep],fabs(m_bar[Tstep]),Chi[Tstep]);
        
        printf("%.2lf%%\tT=%.3lf\t<H>=%.2lf\tC=%.4lf\tm=%.4lf\tChi=%.4lf\n",
               (double)(Tstep) / (double)(Tm) * 100.0, T, H_bar[Tstep], C[Tstep], fabs(m_bar[Tstep]), Chi[Tstep]);
    }
    fclose(fp);
}
