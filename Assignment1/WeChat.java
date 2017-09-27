import java.util.Scanner;

public class WeChat {
	public static void main(String[] args) {
		Scanner input = new Scanner(System.in);
		double total;//红包总金额
		int num;//红包个数
		while(true) {// 输入满足要求的红包总金额（total）及数目（numb）
			System.out.print("请输入红包金额(保留两位小数)：");
			total = input.nextDouble();
			System.out.print("请输入红包份数：");
			num =input.nextInt();
			if(total<=0) {
				System.out.println("请输入大于零的红包金额");
			}
			else if(num<=0){
				System.out.println("请输入大于零的红包份数");
			}
			else if(num>total*100) {
				System.out.println("红包份数过多，请重新来过");
			}
			else {
				break;
			}
		}	
		
		//产生红包金额的方式为，对于每个红包赋予权重（产生一个随机数）
		//然后该红包金额占总的比例为对应随机数占随机数总和的比例
		
		int i;
		double moneynow;//当前分配金额
		int lucky=1;//运气王
		double luckymoney=0.0;//运气王的金额
		double[] ran=new double[num];//各红包的权重
		double sum=0.0;//总权
		double[] money=new double[num];//每个红包内的钱数
		money[num-1]=total;//最初最后一个红包剩的值
		for(i=0;i<num;i++) {//产生每个红包内金额的权重
			ran[i]=Math.random()+0.0001;//防止sum为0
			sum+=ran[i];
		}
		for(i=0;i<num-1;i++) {//产生红包的值
			money[i]=(int)(100.0*ran[i]/sum*(total-0.01*(num-1)))/100.0+0.01;//保证产生两位小数且不为0
			money[num-1]-=money[i];
		}
		
		i=(int)(Math.random()*num);//因最后一个的运气王概率和前面不等，所以随机调换
		sum=money[num-1];
		money[num-1]=money[i];
		money[i]=sum;
		
		for(i=0;i<num;i++) {//输出运气王的情况
			moneynow=money[i];
			System.out.println("第"+ (i+1)+ "份红包的金额为"+ String.format("%.2f", moneynow));//输出该红包的值
			if(moneynow>luckymoney) {//更新运气王的情况
				lucky=i+1;
				luckymoney=moneynow;
			}
		}
		System.out.println("\n 运气王为第"+ lucky + "位，其金额为"+String.format("%.2f", luckymoney));
	}
}
