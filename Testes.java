import java.util.ArrayList;
import java.util.Arrays;


public class Testes{
        
    public static void println(String s){
        System.out.println(s);
    }

	public static void main(String[] args) {
        double[][] v1 = {{1},{2},{3}};
        double[][] v2 = {{0.5},{0.1},{1}};
        double[][] v1_T =LinAlg.transpose(v1); 
        println("V1: ");
        LinAlg.print(v1);
        println("V2");
        LinAlg.print(v1_T);

        double[][] A ={{2,2,2,13},{2,1,0.7,0.57},{11,5,9,2.3},{27,31,1,4}};
        double[][] A_ = LinAlg.inv(A);

        
        double[][] h = LinAlg.hadamard(v1,v2);
        double dot = LinAlg.dot(v1,v2);
        double[][] dotM = LinAlg.matmul(v1_T, v2);
        double[][] I3 = LinAlg.eye(3);
        I3 = LinAlg.prod(I3,3);
        double [][]B = {{3,3,3}};
        LinAlg.print(B);

        LinAlg.print(A);
        LinAlg.print(A_);
        LinAlg.print(LinAlg.matmul(A,A_));

        LinAlg.print(LinAlg.upHessenberg(A));
        LinAlg.QR(A);
        double [][]Q = LinAlg.getQ();
        double [][]R = LinAlg.getR();
        LinAlg.print(Q);
        LinAlg.print(R);

        LinAlg.print(LinAlg.matmul(Q,R));
		
	}

};