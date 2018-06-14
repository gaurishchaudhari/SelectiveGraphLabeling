package iitb.shared.optimization;

import iitb.shared.ArrayUtils;
import iitb.shared.RobustMath;

import java.util.Arrays;

/**
 * finds min_{x \in Delta^n} f(x) using the mirror descent algorithm.
 */

public class MirrorDescentTrainer implements ConvexOptimization {
	
	private int debugLvl;
	private int maxIters;
	double grad[];
	int icall;
	private double epsForConvergence;
	
	public MirrorDescentTrainer(int debugLvl, int maxIter, double epsForConvergence) {
		this.debugLvl = debugLvl;
		this.maxIters=maxIter;
		this.epsForConvergence = epsForConvergence;
	}
	
	public double optimize(double x[], GradientCompute gc) throws Exception {
		Arrays.fill(x,1.0/x.length);
		grad = new double[x.length];
		double p[] = new double[x.length];
		double obj = 0;
		icall = 0;
		double delta = 1e-16;
		double sqrtLnN = Math.sqrt(Math.log(x.length));
		double prevObj = Double.MAX_VALUE;
		do {
			obj = gc.computeFunctionGradient(x, grad);
			double maxNorm = ArrayUtils.normInfty(grad);
			if(debugLvl > 0) {
				System.out.println("Iter "+icall + " obj "+obj + " gradMaxNorm " + maxNorm);
			}
			if (obj < prevObj && prevObj - obj < epsForConvergence*prevObj)
				break;
			if (maxNorm < epsForConvergence) break;
			double stepSize = sqrtLnN/maxNorm/Math.sqrt(icall+1);
			if (obj > prevObj) {
				stepSize = stepSize/2;
			}
			prevObj = obj;
			for (int i = 0; i < p.length; i++) {
				p[i] = 1 + Math.log(x[i] + delta) - stepSize*grad[i];
			}
			double norm = RobustMath.logSumExp(p);
			for (int i = 0; i < p.length; i++) {
				x[i] = Math.exp(p[i]-norm);
				assert(!Double.isNaN(x[i]) || !Double.isInfinite(x[i]));
			}
			
		} while ((icall++ <= maxIters));
		return obj;
	}
	
	public double[] lastGrad() {
		return grad;
	}
	
	@SuppressWarnings("unused")
	public static void main(String args[]) throws Exception {
		final double d[] = new double[]{-2, 3, -1};
		double x[] = new double[d.length];
		MirrorDescentTrainer md = new MirrorDescentTrainer(1, 100, 1e-12);
		double obj = md.optimize(x, new GradientCompute() {
			@Override
			public double computeFunctionGradient(double[] lambda, double[] grad)
					throws Exception {
				double obj = 0;
				for (int i = 0; i < lambda.length; i++) {
					grad[i] = -d[i]*d[i]/lambda[i]/lambda[i];
					obj += d[i]*d[i]/lambda[i];
				}
				return obj;
			}
		});
		System.out.println(Arrays.toString(x));
	}
}
