#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define N 16 //N个粒子
#define LL 4 //N个粒子
int msteps = 1000; //达到热平衡后的模拟次数
int istep =10000;
int nbins=1;
int theta[N+1]; //自旋
double m0; //总自旋
int n;//算符个数
int nb=2*N;
int L;//截断阶数 
//int p,l;//第p行第l个
int v;//v=4*p+l 
int v_first[N+1],v_last[N+1];
int bsites[2][2*N+1];//键的两端 
double beta=32.0; //T (J / kB)wa 
int x[4*80];
int s[80]; 
//int *x=(int*)malloc(4*400*sizeof(int));//所连接的节点 
//int *s=(int*)malloc(400*sizeof(int));//s=2*a+b-1	
double enrg1=0.0,enrg2=0.0,amag1=0.0,amag2=0.0,asusc=0.0,stiff=0.0,ususc=0.0;
double data1[7],data2[7]; 
double wdata1[7],wdata2[7];
//J = 1 

//初始化 
void init(){
	int i,j,q,p;
    L=80;
    for (p = 1; p <= N; p++){
        theta[p] = (rand() % 2) ? 1 : -1;
    }
    for (p = 0; p < L; p++){
        s[p]=0;
    }
    n=0;
    for(j=0;j<LL;j++) {
    	for(i=0;i<LL;i++){
    		q=1+i+j*LL;
    		bsites[0][q]=q;
    		bsites[1][q]=1+(i+1)%LL+j*LL;
    		bsites[0][q+N]=q;
    		bsites[1][q+N]=1+i+LL*((j+1)%LL);
		}
	} 
}

//改变(i,j)处自旋方向能量的改变
void diagonalupdate(){
	int b,op,p;
	for(p=0;p<L;p++){
		op=s[p]; 
		if(op==0){
			b=(rand())%(2*N)+1;
			if(theta[bsites[0][b]]!=theta[bsites[1][b]]){
				if((double)(rand()) / (double)(RAND_MAX)*(double)(L-n)<(double)(N)*beta){
					s[p]=2*b;
					n=n+1;
				}
			}
		}
		else if(op%2==0){
			if(((double)(rand()) / (double)(RAND_MAX)*(double)(N)*beta<(double)(L-n+1))){
				s[p]=0;
				n=n-1;
			}
		}
		else{
			b=op/2;
			theta[bsites[0][b]]=-theta[bsites[0][b]];
			theta[bsites[1][b]]=-theta[bsites[1][b]];
		}
	} 
} 

int flipbit(int ms){
	int ss;
	if(ms%2==0){
		ss=ms+1;
	}
	else{
		ss=ms-1;	
	}	
	return ss;
}

void loopupdate(){
	int p,k,v0,b,i1,i2,v1,v2,f,op;
	for(p=0;p<=N;p++){
		v_first[p]=-1;
		v_last[p]=-1;
	} 
	for(v0=0;v0<4*L;v0=v0+4){
		op=s[v0/4];
		if(op!=0){
			b=op/2;
			i1=bsites[0][b];
			i2=bsites[1][b];
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
			x[v0]=0;
			x[v0+1]=0;
			x[v0+2]=0;
			x[v0+3]=0;
		}
	}

	for(k=1;k<=N;k++){
		f=v_first[k];
		if(f!=-1){
			v2=v_last[k];
			x[v2]=f;
			x[f]=v2;
		}	
	}
	
	for(v0=0;v0<4*L;v0=v0+2){
		if(x[v0]<1){
			continue;
		}
		v1=v0;
		if((double)(rand()) / (double)(RAND_MAX)<0.5){
			while(1){ 
          		s[v1/4]=flipbit(s[v1/4]);
          		x[v1]=-1;
          		v2=flipbit(v1);
				v1=x[v2];
				x[v2]=-1;
          		if (v1==v0) 
				  break;	
			}
		}
		else{
			while(1){ 
          		x[v1]=0;
          		v2=flipbit(v1);
        		v1=x[v2];
        		x[v2]=0;
          		if (v1==v0) 
				  break;	
			}	
		}
		for(k=1;k<=N;k++){
			if(v_first[k]!=-1){
				if(x[v_first[k]]==-1){
					theta[k]=-theta[k];
				} 
			}
			else if((double)(rand()) / (double)(RAND_MAX)<0.5){
					theta[k]=-theta[k];
				}
		}
	}
}

//void adjustcutoff(int step){
//	int mmnew=n+n/3;
//	int kkk;
//	if (mmnew>L){
//		s=(int*)realloc(s,mmnew*sizeof(int)); 
//		for(kkk=L;kkk<mmnew;kkk++)
//			s[kkk]=0;
//		L=mmnew;
//		x=(int*)realloc(x,4*L*sizeof(int)); 
//	}
//}
//
//double mmm000(int theta[N]){
//	double mmmm=0.0; 
//	int ka;
//	for(ka=0;ka<N;ka++){
//		if((ka%LL+ka/LL)%2==0){
//			mmmm+=(double)(theta[ka]);
//		}
//		else{
//			mmmm-=(double)(theta[ka]);
//		}
//	}
//	mmmm/=2.0;
//	return mmmm;
//}

int main(){
	FILE *fp;
	int i,p;
 	char name[30];
 	printf("请输入文件名  :");
 	scanf("%s",name);
	fp = fopen(name, "w+");
	int Tstep,k,kq,step,b,kk,sum;
	double am1,am2,ax1;
	double jj[2];
    for(i=0;i<7;i++){
		data1[i]=0.0;
		data2[i]=0.0;
	}
	srand((int)time(0));
	enrg1=0.0;enrg2=0.0;amag1=0.0;amag2=0.0;asusc=0.0;stiff=0.0;ususc=0.0;
    init();
	for(k=1;k<=istep;k++){
		diagonalupdate();
		loopupdate();
	//	adjustcutoff(k); 	 
	}
            //热平衡
            //模拟ustep次
    for(kq=0;kq<nbins;kq++){
		for (step = 0; step < msteps; step++){
			diagonalupdate();
			loopupdate(); 	 		
			m0=0.0; 
			for(k=1;k<=N;k++){
				if(((k-1)%LL+(k-1)/LL)%2==0){
					m0+=(double)(theta[k]);
				}
				else{
					m0-=(double)(theta[k]);
				}
			}
			m0/=2.0;
//			
//			printf("n=%d\tL=%d\n",n,L);
//			for(k=0;k<L;k++){
//				printf("s[%d]=%d\n",k,s[k]);
//			} 
//			
			am1=0.0;
			am2=0.0;
			ax1=0.0;
			jj[0]=0.0;
			jj[1]=0.0;
			if(step==5){
						for(k=1;k<=N;k++){
							printf("theta[%d]=%d| ",k,theta[k]);	
						}
						printf("\n  m0=%lf\n",bsites[0][b],bsites[1][b],m0);
					}
			for(p=0;p<L;p++){
				if(s[p]==0){
					continue;
				}
				else if(s[p]%2==1){
					b=s[p]/2;
					if(step==5){
						for(k=1;k<=N;k++){
							printf("theta[%d]=%d| ",k,theta[k]);	
						}
						printf("\n b1=%d b2=%d | m0=%lf\n",bsites[0][b],bsites[1][b],m0);
					}
//					printf("b is %d,theta is %d %d, m0 is %lf\n",b,theta[bsites[0][b]],theta[bsites[1][b]],m0);
					theta[bsites[0][b]]=-theta[bsites[0][b]];
					theta[bsites[1][b]]=-theta[bsites[1][b]];
					jj[(b-1)/N]+=(double)(theta[bsites[1][b]]);
//					m0=mmm000(theta);
					if(((bsites[0][b]-1)%LL+(bsites[0][b]-1)/LL)%2==0){
						m0+=2.0*(double)(theta[bsites[0][b]]);
					}
					else{
						m0-=2.0*(double)(theta[bsites[0][b]]);
					}
				}
				ax1+=m0;
				am1+=fabs(m0);
				am2+=m0*m0;
			}
			if(n!=0){
				ax1=(ax1*ax1+am2)/((double)(n)*(double)(n+1));
				am1/=(double)(n);
				am2/=(double)(n);	
			}
			else{
				am1=fabs(m0);
				am2=m0*m0;
				ax1=am2;
			}
//			printf("am1是%lf,n是%d\n",am1,n);
//			return 1;
			enrg1+=(double)(n);
 			enrg2+=(double)(n)*(double)(n);
 			amag1+=am1;
 			amag2+=am2;
			asusc+=ax1;
		 	stiff+=0.5*(jj[0]*jj[0]+jj[1]*jj[1]);
		 	sum=0;
			for(kk=1;kk<=N;kk++){
		 		sum+=theta[kk];
			}
			ususc+=(double)(sum/2)*(double)(sum/2);
        }	
        enrg1/=msteps;
 		enrg2/=msteps;
 		amag1/=msteps;
 		amag2/=msteps;
 		asusc/=msteps;
 		stiff/=msteps;
 		ususc/=msteps;
 		enrg2=(enrg2-enrg1*(enrg1+1.0))/(double)(N);
 		enrg1=enrg1/beta/(double)(N)-0.5;
 		amag1/=(double)(N);
 		amag2/=(double)(N);
 		asusc/=(double)(N)/beta;    
 		ususc/=(double)(N)/beta;    
 		stiff/=(double)(N)*beta;
	//	if(kq%2000==1){
	//		printf("%lf\n\n",enrg1);
	//	}
 		data1[0]+=enrg1;
 		data1[1]+=enrg2;
 		data1[2]+=amag1;
 		data1[3]+=amag2;
 		data1[4]+=asusc;
 		data1[5]+=stiff;
 		data1[6]+=ususc;
	//	if(kq%2000==1){
	//		printf("%lf\n\n",data1[0]);
	//	}
 		data2[0]+=enrg1*enrg1;
 		data2[1]+=enrg2*enrg2;
 		data2[2]+=amag1*amag1;
 		data2[3]+=amag2*amag2;
 		data2[4]+=asusc*asusc;
 		data2[5]+=stiff*stiff;
 		data2[6]+=ususc*ususc;
		for(i=0;i<7;i++){
			wdata1[i]=data1[i]/(double)(kq+1);
   			wdata2[i]=data2[i]/(double)(kq+1);
    		wdata2[i]=sqrt(fabs(wdata2[i]-wdata1[i]*wdata1[i])/(double)(kq+1));
		}
	//	if(kq%100==1){
			printf("-E/N=%.4lf\t%.4lf\tC/N=%.4lf\t%.4lf\t<|m|>=%.4lf\t%.4lf\tS(pi,pi)=%.4lf\t%.4lf\tX(pi,pi)=%.4lf\t%.4lf\trho_s=%.4lf\t%.4lf\t\tX(0,0)=%.4lf\t%.4lf\n",wdata1[0],wdata2[0],wdata1[1],wdata2[1],wdata1[2],wdata2[2],wdata1[3],wdata2[3],wdata1[4],wdata2[4],wdata1[5],wdata2[5],wdata1[6],wdata2[6]);
      		fprintf(fp,"-E/N=%.4lf\t%.4lf\tC/N=%.4lf\t%.4lf\t<|m|>=%.4lf\t%.4lf\tS(pi,pi)=%.4lf\t%.4lf\tX(pi,pi)=%.4lf\t%.4lf\trho_s=%.4lf\t%.4lf\t\tX(0,0)=%.4lf\t%.4lf\n",wdata1[0],wdata2[0],wdata1[1],wdata2[1],wdata1[2],wdata2[2],wdata1[3],wdata2[3],wdata1[4],wdata2[4],wdata1[5],wdata2[5],wdata1[6],wdata2[6]);
	//	}
		enrg1=0.0;
		enrg2=0.0;
		amag1=0.0;
		amag2=0.0;
		asusc=0.0;
		stiff=0.0;
		ususc=0.0;
	}
	fclose(fp);
}
