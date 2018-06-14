package iitb.shared.optimization;

public interface GradientCompute {
    double computeFunctionGradient(double lambda[], double grad[])  throws Exception ;
}
