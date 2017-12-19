#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N 16 //N������
int Tm = 100; //T�����
double dT = 0.05; //T����
int ustep = 10000; //�ﵽ��ƽ����ģ�����
int mstep =1000;
double eps = 1e-8;
int theta[N]; //����
int m0; //������
int n;//�������
int L=N/3+N;//�ضϽ��� 
int p,l;//��p�е�l��
int v;//v=4*p+l 
int x[4*(N/3+N)];//�����ӵĽڵ� 
int s[N/3+N];//s=2*a+b-1
int v_first[N],v_last[N];
int bsites[2][2*N];//�������� 
double T = 6.0; //T (J / kB)wa 
//J = 1 

//��ʼ�� 
void init(){
	int i,j,q;
    m0 = 0;
    for (p = 0; p < L; p++){
        theta[p] = (rand() % 2) ? 1 : -1;
        s[p]=0;
    }
    n=0;
//    for (p = 0; p <4*L; p++){
//    	x[p]=-1;
//    }
    for(j=0;j<4;j++) {
    	for(i=0;i<4;i++){
    		q=i+j*4;
    		bsites[0][q]=q;
    		bsites[1][q]=(i+1)%4+j*4;
    		bsites[0][q+16]=q;
    		bsites[1][q+16]=i+4*((j+1)%4);
		}
	} 
}

//�ı�(i,j)���������������ĸı�
void diagonalupdate(){
	int b;
	for(p=0;p<L;p++){
		if(s[p]==0){
			b=rand() % N+1;
			if(theta[bsites[0][b-1]]==theta[bsites[1][b-1]]){
				continue;
			}
			if((double)(rand()) / (double)(RAND_MAX)*(double)(L-n)<(double)(N)/(double)(T)/2){
				s[p]=2*b;
				n=n+1;
			}
		}
		else if(s[p]%2==0){
			if(((double)(rand()) / (double)(RAND_MAX)*(double)(N)/2<(double)(L-n+1)*(double)(T))){
				s[p]=0;
				n=n-1;
			}
		}
		else{
			b=s[p]/2;
			theta[bsites[0][b-1]]=-theta[bsites[0][b-1]];
			theta[bsites[1][b-1]]=-theta[bsites[1][b-1]];
		}
	} 
} 

int flipbit(int s){
	int ss;
	if(s%2==0){
		ss=s+1;
	}
	else{
		ss=s-1;	
	}	
	return ss;
}

void loopupdate(){
	int k,v0,b,i1,i2,v1,v2,f;
	for(p=0;p<L;p++){
		v_first[p]=-1;
		v_last[p]=-1;
	} 
	for(v0=0;v0<4*L;v0=v0+4){
		p=v0/4;
		if(s[p]!=0){
			b=s[p]/2;
			i1=bsites[0][b-1];
			i2=bsites[1][b-1];
			v1=v_last[i1];
			v2=v_last[i2];
			if(v1!=-1){
				x[v1]=v0;
				x[v0]=v1;
			}
			else{
				v_first[i1]=v0;
			}
			if(v2!=-1){
				x[v2]=v0+1;
				x[v0+1]=v2;
			}
			else{
				v_first[i2]=v0+1;
			}
			v_last[i1]=v0+2;
			v_last[i2]=v0+3;
		}
		else{
			x[v0]=-1;
			x[v0+1]=-1;
			x[v0+2]=-1;
			x[v0+3]=-1;
		}
	}

	for(k=0;k<N;k++){
		f=v_first[k];
		if(f!=-1){
			v2=v_last[k];
			x[f]=v2;
			x[v2]=f;
		}	
	}
	
	for(v0=0;v0<4* L;v0=v0+2){
		if(x[v0]<0){
			continue;
		}
		v1=v0;
		if((double)(rand()) / (double)(RAND_MAX)<0.5){
			while(1){ 
          		s[v1/4]=flipbit(s[v1/4]);
          		x[v1]=-2;
          		v2=flipbit(v1);
				v1=x[v2];
				x[v2]=-2;
          		if (v1==v0) 
				  break;	
			}
		}
		else{
			while(1){ 
          		x[v1]=-1;
          		v2=flipbit(v1);
        		v1=x[v2];
        		x[v2]=-1;
          		if (v1==v0) 
				  break;	
			}	
		}
		for(k=0;k<N;k++){
			v=v_first[k]; 
			if(v==-1){
				if((double)(rand()) / (double)(RAND_MAX)<0.5){
					theta[k]=-theta[k];
				}
				else if(x[v]==-2){
					theta[k]=-theta[k];
				}
			}
		}
	}
}

int main(){
    double m_square_bar[Tm]; //�žغʹž�ƽ���ľ�ֵ

    FILE *fp;
		int i;
 		char name[30];
 		printf("�������ļ���  :");
 		scanf("%s",name);
		fp = fopen(name, "w+");
	int Tstep,k,step,b;
	double t[250];
	double tau=0,kk=0,fang=0;
	double m2,mm;
	double m_bar[ustep];
    //for (Tstep = 1; Tstep <= Tm; Tstep++, T += dT){
        m_square_bar[Tstep] = 0;

        init();
		for(k=0;k<mstep;k++){
			diagonalupdate();
			loopupdate();
		}
            //��ƽ��
            //ģ��ustep��
            for (step = 0; step < ustep; step++){
                for(k=0;k<5000;k++){
                	diagonalupdate();
					loopupdate();
				}
				m0=0; 
				for(k=0;k<N;k++){
					m0+=(-1)^(k%4+k/4)*theta[k];
				} 
				m2=0;
				for(p=0;p<L;p++){
					if(s[p]%2==1){
						b=s[p]/2;
						theta[bsites[0][b-1]]=-theta[bsites[0][b-1]];
						theta[bsites[1][b-1]]=-theta[bsites[1][b-1]];
						m0+=2*(-1)^((bsites[0][b-1])%4+(bsites[0][b-1])/4)*theta[bsites[0][b-1]];
					}
					if(s[p]!=0)
						m2+=(double)(m0*m0);
				}
				if(n==0){
					m_square_bar[Tstep] += m2/(double)(N*N);	
				}
				else{
					m_square_bar[Tstep] += m2/(double)(n*N*N);
				}
				m_bar[step] =fabs((double)(m0));
                k += fabs((double)(m0));
                fang += m0*m0;
            }
			k /= (double)(ustep) * (double)(N);
            fang /=(double)(ustep) * (double)(N * N);
            for(i=1;i<251;i++){
            	tau=0;
            	for (step = 0; step < ustep; step++){
            		tau+=m_bar[step]*m_bar[(step+i)%ustep];
            	}
            	tau /= (double)(ustep) * (double)(N * N);
							t[i-1]=(tau-k*k)/(fang-k*k);
            
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
