package iitb.shared.optimization;

import iitb.shared.ArrayUtils;
import riso.numerical.LBFGS;

public class LBFGSTrainer implements ConvexOptimization {
	
	public int icall;
	double grad[];
	int debugLvl, mForHessian, maxIters;
	double epsForConvergence;
	public boolean minimize=false;
	double f;
	
	public LBFGSTrainer(int debugLvl, int mForHessian, int maxIter, double epsForConvergence) {
		this(debugLvl,mForHessian,maxIter,epsForConvergence,false);
	}
	
	public LBFGSTrainer(int debugLvl, int mForHessian, int maxIter, double epsForConvergence, boolean minimize) {
		this.debugLvl = debugLvl;
		this.mForHessian = mForHessian;
		this.maxIters=maxIter;
		this.epsForConvergence = epsForConvergence;
		this.minimize=minimize;
	}
	
	public double optimize(double lambda[], GradientCompute gc) throws Exception {
		return doTrain(lambda, gc);
	}
	
	public double doTrain(double lambda[], GradientCompute gc) throws Exception {
		try {
			doTrainE(lambda, gc);
		}
		catch (LBFGS.ExceptionWithIflag e)  {
			System.err.println( "lbfgs failed.\n"+e );
			if (e.iflag == -1) {
				System.err.println("Possible reasons could be: \n \t 1. Bug in the feature generation or data handling code\n\t 2. Not enough features to make observed feature value==expected value\n");
			}
			return f;
		}
		return f;
	}
	
	public double doTrainE(double lambda[], GradientCompute gc) throws Exception,  LBFGS.ExceptionWithIflag{
		f=0;
		double xtol = 1.0e-16; // machine precision
		int iprint[] = new int [2], iflag[] = new int[1];
		icall=0;
		grad = new double[lambda.length];
		double diag[] = new double[lambda.length];
		iprint [0] = debugLvl-2;
		iprint [1] = debugLvl-1;
		iflag[0]=0;
		do {
			f = gc.computeFunctionGradient(lambda,grad);
			if (!minimize) {
				f = -1*f; // since the routine below minimizes and we want to maximize logli
				for (int j = 0 ; j < lambda.length ; j ++) {
					grad[j] *= -1;
				}
			}
			if (debugLvl > 0) {
				System.out.println("Iter "+icall + " f "+f  + " gradNorm "+ ArrayUtils.norm2Square(grad) + " xnorm "+ArrayUtils.norm2Square(lambda)); // + " x: "+Arrays.toString(lambda)
			}
			
			LBFGS.lbfgs(lambda.length, mForHessian, lambda, f, grad, false, diag, iprint, 
					epsForConvergence, xtol, iflag);
			icall += 1;
			
		} while (( iflag[0] != 0) && (icall <= maxIters));
		return f;
	}
	
	public double[] lastGrad() {
		return grad;
	}
	
	public double lastObjective() {
		return f;
	}
}