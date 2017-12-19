#include<stdio.h>









long long int f(int a,int b){
	long long int c;
	if(b==1){
		return a;
	}
	if(b%2==1){
		return a*f(a,b-1);
	}
	else{
		c=f(a,b/2);
		return c*c;
	}
}


int main(){
	int a,b;
	long long int c; 
	scanf("%d %d",&a,&b);
	c=f(a,b);
	printf("%lld",c);	
	 
}


