import java.util.ArrayList;
import java.lang.Math;
import java.util.Random;
import java.util.Arrays;
import java.lang.System;

/*
Autor: Vitor Barbosa
Implementação baseada nos algoritmos disponíveis em:
Blackledge, J.: Digital Signal Processing (Second Edition). Horwood Publishing, vol: ISBN: 1-904275-26-5. 2006.
*/

public class LinAlg{

	//Variáveis temporárias usadas para retornar múltiplos valores. Não é um modo ótimo, mas funciona
	private static double[][] _Q;
	private  static double[][] _R;

	public static double[][] getQ(){return copy(_Q);}
	public static double[][] getR(){return copy(_R);}

	public static void print(int[] A) {System.out.println(Arrays.toString(A));}
	public static void print(double A) {System.out.println(A);}
	public static void print(double[] A){System.out.println(Arrays.toString(A));}

	public static void print(double[][] A){
		String s = "";
		for(double [] row:A) s+= " \n"+Arrays.toString(row);
		s+="\n";
		System.out.println(s);
	}

	public static int[] dim(double[][] A){
		int r = A.length;
			if(r==0) return new int[]{-1,-1};
		int c = A[0].length;
		return new int[]{r,c};
	}

	public static double[][] transpose(double [][] A){
		int c = A[0].length, r = A.length;
		double[][] AT = new double[c][r];
		for(int i=0;i<r;i++) for(int j=0;j<c;j++) AT[j][i] = A[i][j];
		return AT;
	}

	public static double[][] copy(double [][] A){
		int c = A[0].length, r = A.length;
		double[][] B = new double[r][c];
		for(int i=0;i<r;i++) for(int j=0;j<c;j++) B[i][j] = A[i][j];
		return B;
	}

	public static double[] copy(double[] A){
		double [] vec = new double[A.length];
		for(int i=0;i<A.length;i++) vec[i] = A[i];
		return vec;
	}

	public static double dot(double[][] A, double[][] B){
		double sum = 0.0;
		int i=0,j=0;
		int [] dimA = dim(A);
		if(dimA[0]>=dimA[1]) // Vetores coluna
			for(;i<dimA[0];i++) sum+=A[i][0]*B[i][0];
		else 
			for(;i<dimA[1];i++) sum+=A[0][i]*B[0][i];
			//for(;i<dim(A)[0];i++) sum+=A[i][0]*B[0][i];
		return sum;
	}
	//Multiplica matriz por escalar
	public static double [][] prod(double k, double[][]A){return prod(A,k);}
	public static double [][] prod(double[][] A,double k){
		int c = A[0].length, r = A.length;
		double[][] B = new double[r][c];
		for(int i=0;i<r;i++) for(int j=0;j<c;j++) B[i][j] = k*A[i][j];
		return B;
		
	}
	public static double mean(double[]A){
		double sum = 0;
		for(double x: A) sum+=x;
		return sum/A.length;
	}

	public static double stdDev(double[] A){
		double sum = 0, mu = mean(A);
		for(double x:A) sum+=(x-mu)*(x-mu);
		return Math.sqrt(sum/(A.length-0.0));
	}

	//Retorna o conjunto diferença entre os conjuntos A e B
	public static int[] setDifference(int[]A, int []B){
		int n = A.length, i=0;
		for(int a:A) for(int b:B) if(b==a) n--; // Subtrai os elementos comuns aos dois conjuntos da contagem
		int []diff = new int[n];
		for(int a:A){
			boolean shared = false;
			for(int b:B) if(b==a) shared = true;
			if(!shared) diff[i++] = a;
		}
		return diff;
	}

	//Multiplica até 5 matrizes. A recursão não é eficiente, mas facilita muito o uso
	public static double[][] matmul(double[][]A, double[][] B, double[][] C, double[][] D, double[][] E){return matmul(matmul(A,B,C,D),E);}
	public static double[][] matmul(double[][]A, double[][] B, double[][] C, double[][] D){return matmul(matmul(A,B,C),D);}
	public static double[][] matmul(double[][]A, double[][] B, double[][] C){return matmul(matmul(A,B),C);}
	public static double[][] matmul(double[][] A, double[][] B){
		int rowsA = A.length, rowsB = B.length;
		int colsA = A[0].length, colsB = B[0].length;
		//Dimensões: A(mxp) * B(pxn) = X(mxn)
		double[][] X = new double[rowsA][colsB];
		double sum=0;
		for(int i=0;i<rowsA;i++){
			int j=0;
			for(;j<colsB;j++){
				sum=0;
				for(int k=0;k<colsA;k++)
					sum+=A[i][k]*B[k][j];
				X[i][j] = sum;
			}
		}
		return X;
	}
	
	//Retorna a matriz identidade
	public static double[][] eye(int dim){
		double[][] I = new double[dim][dim];
		for(int i =0;i<dim;i++)
			for(int j = 0;j<dim;j++)
				if(i==j) I[i][i] = 1;
		return I;
	}

	public static double[][] zeros(int r, int c){return new double[r][c];}
	public static double[][] zeros(int r){return new double[r][1];}

	//Retorna array de int com valores de uma PA
	public static int[] linspaceInt(int start, int end, int step){
		int [] seq = new int[(end-start)/step];
		seq[0] = start;
		for(int n=0;n<end;n++) seq[n] = seq[0] + n*step;
		return seq;
	}
	//Retorna array de double valores de uma PA
	public static double[] linspace(int start, int end, int step){
		double [] seq = new double[(end-start)/step];
		seq[0] = start;
		for(int n=0;n<end;n++) seq[n] = seq[0] + n*step;
		return seq;
	}

	//Retorna true se a array contém o elemento
	public static boolean contains(int [] A, int x){
		for(int el:A) if(x==el) return true;
		return false;
	}
	

	public static double[] rowVec (double[][] A, int row){
		//public static void arraycopy(Object src, int srcPos, Object dest, int destPos, int length)
		double[] V = new double[A[0].length];
		System.arraycopy(A[row],0,V,0,A[row].length);
		return V;
	}

	//Extrai vetor coluna da matriz
	public static double[] colVec(double[][]A, int col){
		double[] C = new double[A.length];
		for(int v=0;v<A.length;v++) C[v] = A[v][col];
		return C;
	}

	//Extrai vetor coluna em formato 2d da matriz
	public static double[][] colMat(double[][]A, int col){
		double[][] C = new double[A.length][1];
		for(int v=0;v<A.length;v++) C[v][0] = A[v][col];
		return C;
	}


	public static void copyCol(double[][]src, int srcCol, double[][]dst, int dstCol){
			for(int i=0;i<src.length;i++) dst[i][dstCol] = src[i][srcCol];
	}

	//Retorna matriz diagonal com os elementos do vetor d
	public static double[][] diag(double[] d){
		double[][] D = new double[d.length][d.length];
		for(int z=0;z<d.length;z++) D[z][z] = d[z];
		return D;
	}

	//Gera um histograma com 'bins' divisoes
	public static double[] hist(double[] x, int bins){
		double max=x[0], min=x[0];
		double []h = new double[bins];
		//Pega maior e menor elementos
		for(double el: x)	{if(el>max)max=el; else if(el<min) min=el;}
		//Normaliza e incrementa os bins apropriados
		for(double el: x){
			double div = (max-min)/bins;
			int n = (int)((el-min)/div);
			n = Math.min(Math.max(0,n),bins-1);
			h[n]++;
		} 
		return h;
	}


	//Computa a inversa pelo algoritmo de Jordan
	//A matriz precisa ser quadrada
	public static double[][] inv(double[][] Mat){
		int dim = Mat.length;
		double[][] A = copy(Mat);
		double[][] B = eye(dim);
		//Passo k
		for(int k =0;k<dim;k++)
		{
			//Linhas
			for(int i=0;i<dim;i++){
				if(i!=k && A[i][k]!=0){
					double scale = A[i][k]/A[k][k];
					//Colunas
					for(int j=0;j<dim;j++){
						if(j>=k) A[i][j] = A[i][j] - scale*A[k][j];
						B[i][j] = B[i][j] - scale*B[k][j];
					}
				}
			}
		}
		//Se i==k (linha do passo, deve ser 1 0 0 por exemplo)
		for(int i=0;i<dim;i++){
			double scale = 1.0/A[i][i];
			for(int j=0;j<dim;j++){
				A[i][j] = A[i][j]*scale;
				B[i][j] = B[i][j]*scale;
			}
		}
		//print(A);
		return B;
		
	}
			
	public static double [] hadamard(double[] A, double[] B){
		double []X = new double[A.length];
		for(int h=0;h<A.length;h++) X[h] = A[h]*B[h];
		return X;
	}

	//Retorna o produto de hadamard (elemento a elemento)
	public static double[][] hadamard(double[][] A, double[][] B){
		int rows = A.length;
		int cols = A[0].length;
		double[][] X = new double[rows][cols];
		for(int i=0;i<rows; i++)
			for(int j=0;j<cols;j++)
				X[i][j] = A[i][j]*B[i][j];
		return X;
	}

	public static double[][] sub(double[][] A, double[][] B){
			int rows = A.length;
			int cols = A[0].length;
			double[][] X = new double[rows][cols];
			for(int i=0;i<rows; i++)
				for(int j=0;j<cols;j++)
					X[i][j] = A[i][j]-B[i][j];
			return X;
	}

	public static double[][] add(double[][] A, double[][] B){
			int rows = A.length;
			int cols = A[0].length;
			double[][] X = new double[rows][cols];
			for(int i=0;i<rows; i++)
				for(int j=0;j<cols;j++)
					X[i][j] = A[i][j]+B[i][j];
			return X;
	}

	/*
	* Gera números aleatórios de distribuição normal pela Transformada Box-Muller
	*/
	public static double[][] gauss(double mu, double sigma, int N){
		double [][] X = new double[N][1];
		for(int i=0;i<N;i++) X[i][1] = gauss(mu,sigma);
		return X;
	}
	public static double gauss(double mu, double sigma){
		double a=0,b=0,r=2000;
		while(r>=1.0||r<=0.0){
			a = Math.random();
			b = Math.random(); //0<a<1
			a = 2.0*a-1.0; //-1<a<1
			b= 2.0*b-1.0;
			r = a*a + b*b; //0<r<2
		}
		//r = r/2.0; //0<r<1
		double scale = Math.sqrt(-2.0*Math.log(r)/r);
		double g1 = a*scale; //sigma =1, mu = 0
		double g2 = b*scale;
		//g1 = (g1+g2)/2;
		return mu + sigma*g1; //Converte para mu e sigma desejados.
	}
	
	//Computa a fatoração QR usando rotações de Given
	//Altera as matrizes _Q e _R, que podem ser lidas com getQ() e getR()
	public static void QR(double[][] Mat){
		int m = Mat.length, n = Mat[0].length;
		double [][] A = copy(Mat);
		//A = upHessenberg(A); // Otimização para autovalores, mas não funcionou para o caso geral
		print(A);
		double [][] Q = eye(m);
		double [][] R = new double[m][n];
		double [][] P;
		double theta = 0.0;
		//Percorremos a diagonal de P até o penúltimo elemento
		//A restrição k<n-1 foi adicionada para aceitar matrizes retangulares
		for(int k=0;(k<m-1)&&(k<n-1);k++){ 
			for(int i=k+1;i<m;i++){// Expandimos a matriz de rotação para zerar cada linha sub-diagonal de A
				//System.out.printf("i= %d, k= %d \n",i,k);
				if(A[i][k]!=0){
					theta = -Math.atan(A[i][k]/A[k][k]); //Angulo de rotacao
					//Mat de rotação da forma [c s][-s c] ->[c 0 s][0 1 0][-s 0 c] ...
					P = eye(m);
					P[k][k] = Math.cos(theta);
					P[k][i] = Math.sin(theta);
					P[i][k] = -Math.sin(theta);
					P[i][i] = Math.cos(theta);
					//System.out.println("Passo "+i);
					//print(P);
					A = matmul(transpose(P),A); // (P1*P2*P3...)^T*A1 =(...*P3^T*P2^T*P1^T)*A1 = R 
					Q = matmul(Q,P); //P1*P2*P3*... = Q
				}
			}
		}
		R = A;
		_Q = Q;
		_R = R;
	}

	//Converte uma matriz para a forma de Hessenberg superior, ie [ [1 x x x],[x 1 x x],[0 x 1 x],[0 0 x 1] ] 
	public static double [][]upHessenberg (double [][]Mat ){
		int m = Mat.length;
		double [][]A = copy(Mat);
		double [][]M; 
		double [][]M_;
		for(int k=0;k<m-2;k++){
			M = eye(m);
			M_= eye(m);
			for(int i=k+2;i<m;i++) {
				M[i][k+1] = A[i][k]/A[k+1][k];
				M_[i][k+1] = -A[i][k]/A[k+1][k]; //Inversa
			}

			//A(k+1) = M^-1*A*M
			A = matmul(M_,A); 
			A = matmul(A,M);
			print(A);
		}
		return A;
	}
}


