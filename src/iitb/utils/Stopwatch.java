package iitb.utils;

public class Stopwatch { 

	private long start;

	public Stopwatch() {
		start = System.nanoTime();
	} 
	
	public void restart(){
		start = System.nanoTime();
	}
	
	public double elapsedTime() {
		long now = System.nanoTime();
		return (now - start) / 1000000000.0 ;
	}
	
	public String elapsedTimes() {
		long now = System.nanoTime();
		return ""+(now - start) / 1000000000.0 +"s";
	}
}
