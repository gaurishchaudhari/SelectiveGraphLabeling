package iitb.shared.optimization;

public interface ConvexOptimization {
	public double optimize(double x[], GradientCompute gc) throws Exception;
}
