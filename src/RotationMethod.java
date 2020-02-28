import java.util.ArrayList;
import java.util.List;

public class RotationMethod {
    private double copyMatrA[][];
    private double matrA[][];
    private double matrT[][];
    private double matrU[][];
    private double sinFi;
    private double cosFi;
    private final double epsilon = 0.00001;
    private int n;
    private int count = 0;
    private double polinom[] = {1, -3.1966884499998454, 3.796847573496825, -2.0678062361747767, 0.5082483413018405, -0.044096040836169914};
    private double lambda[];
    private List<double[]> discrepancy;

    RotationMethod(double matrA[][]) {
        n = matrA.length;
        this.copyMatrA = new double[n][n];
        this.lambda = new double[n];
        this.matrT = new double[n][n];
        this.matrA = new double[n][n];
        this.matrU = new double[n][n];
        for(int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                this.matrA[i][j] = 0;
                for (int k = 0; k < n; k++) {
                    this.matrA[i][j] += matrA[k][i] * matrA[k][j];
                }
                this.copyMatrA[i][j] = this.matrA[i][j];
                if(i == j){
                    matrU[i][i] = 1;
                }
                else{
                    matrU[i][j] = 0;
                }
            }
        }
    }
    private Pair maxElement(){
        double max = Math.abs(matrA[0][1]);
        int firstIndex = 0;
        int secondIndex = 1;
        for(int i = 0; i < n; i++){
            for(int j = 1; j < n; j++){
                if(j > i && Math.abs(matrA[i][j]) > max){
                    max = Math.abs(matrA[i][j]);
                    firstIndex = i;
                    secondIndex = j;
                }
            }
        }
        return new Pair(firstIndex, secondIndex);
    }

    private void createNewMatrT(Pair pair) {
        double tg2Fi = 2*matrA[pair.getFirst()][pair.getSecond()] /
                                (matrA[pair.getFirst()][pair.getFirst()] -
                                 matrA[pair.getSecond()][pair.getSecond()]);
        double cos2Fi = 1 / Math.sqrt(1 + Math.pow(tg2Fi, 2));
        cosFi = Math.sqrt((1 + cos2Fi) / 2);
        sinFi = Math.sqrt((1 - cos2Fi) / 2) * Math.signum(tg2Fi);
                                /* matrA[pair.getFirst()][pair.getSecond()] *
                                (matrA[pair.getFirst()][pair.getFirst()] -
                                 matrA[pair.getSecond()][pair.getSecond()]));*/

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    matrT[i][i] = 1;
                } else {
                    matrT[i][j] = 0;
                }
            }
        }
        if (matrA[pair.getSecond()][pair.getFirst()] != 0) {
            matrT[pair.getFirst()][pair.getFirst()] = cosFi;
            matrT[pair.getFirst()][pair.getSecond()] = -sinFi;
            matrT[pair.getSecond()][pair.getSecond()] = cosFi;
            matrT[pair.getSecond()][pair.getFirst()] = sinFi;
        }
    }

    private double difference(){
        double diff = 0;
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(j != i){
                    diff += Math.pow(matrA[i][j], 2);
                }
            }
        }
        return diff;
    }

    private void createNewMatrU(){
        matrU = multiply(matrU, matrT);
    }
    private void createNewMatrA(){
        matrA = multiply(transp(matrT), multiply(matrA, matrT));
    }
    private double[][] transp(double[][] matr){
        double transpMatr[][] = new double[n][n];
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                transpMatr[i][j] = matr[j][i];
            }
        }
        return transpMatr;
    }

    private double[][] multiply(double matrA[][], double matrB[][]){
        double result[][] = new double[n][n];
        for(int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = 0;
                for (int k = 0; k < n; k++)
                    result[i][j] += matrA[i][k] * matrB[k][j];
            }
        }
        return result;
    }


    public void rotationMethod(){
        printPriorityCount();
        while(difference() > epsilon){
            count++;
            createNewMatrT(maxElement());
            createNewMatrA();
            createNewMatrU();
            //printA();
        }
        for(int i = 0; i < n; i++){
            lambda[i] = matrA[i][i];
        }
        System.out.println("count = " + count);
    }



    public void printDiscrepancy(){
        double discrepancy[] = new double[n];
        for(int i = 0; i < n; i++){
            discrepancy[i] = 0;
            for(int j = 0; j < n; j++){
                discrepancy[i] += Math.pow(lambda[i], n-j) * polinom[j];
            }
            discrepancy[i] += polinom[n];
            System.out.println("discrepancy for lambda = " + lambda[i]);
            System.out.format("%25s", discrepancy[i] + "\n");
        }
    }

    private double [] multiply(double[][] matr, double[] vector){
        double result[] = new double[n];
        for(int i = 0; i < n; i++)
        {
            result[i]=0;
            for(int j = 0; j < n; j++)
            {
                result[i] += matr[i][j]*vector[j];
            }
        }
        return result;
    }


    private double[] multiply(double[] vector, double lambda){
        double[] result = new double[n];
        for(int i = 0; i < n; i++){
            result[i] = vector[i] * lambda;
        }
        return result;
    }
    private double[] minus(double[] vectorA, double[] vectorB){
        double[] result = new double[n];
        for(int i = 0; i < n; i++){
            result[i] = vectorA[i] - vectorB[i];
        }
        return result;
    }

    private double[] createVectorDiscrepancy(double[] eigenVector, double lambda){
        return minus(multiply(copyMatrA, eigenVector), multiply(eigenVector, lambda));
    }

    public void createVectorsDiscrepancy(){
        discrepancy = new ArrayList<>();
        for(int i = 0; i < n; i++){
            discrepancy.add(createVectorDiscrepancy(transp(matrU)[i], lambda[i]));
        }
    }

    public void printVectorsDiscrepancy(){
        int count = 0;
        System.out.println("\nEigen vectors discrepancy\n");
        for(double[] vector: discrepancy){
            System.out.println("discrepancy for lambda = " + lambda[count]);
            for(double item : vector){
                System.out.format("%e   ", item);
            }
            System.out.println();
            count++;
        }
    }

    public void printA(){
        System.out.println("Matrix A");
        for (double[] row: matrA){
            for(double item: row){
                System.out.format("%25s", item + "    ");
            }
            System.out.println();
        }
        System.out.println();
    }

    public void printU(){
        System.out.println("Matrix U");
        for (double[] row: matrU){
            for(double item: row){
                System.out.format("%25s", item + "    ");
            }
            System.out.println();
        }
        System.out.println();
    }

    private void printPriorityCount(){
        System.out.println("Priority count = " + (Math.ceil((Math.log(epsilon) - Math.log(difference())) / Math.log(1.0 - 2.0/n/(n-1)))));
    }
}
