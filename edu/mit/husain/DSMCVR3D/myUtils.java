/** (C) 2008 All Rights Reserved 
 * 
 * @author Husain A. Al-Mohssen
 *
 */
package edu.mit.husain.DSMCVR3D;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Random;
import static java.lang.Math.*;

public class myUtils {
	public static Random r=new Random(218232);//random seed :Parameters controls this now
	public static double nextDouble(){return r.nextDouble();}//shorthand
	/** Gaussians with mean 0 and variance=1.0*/
	final public static double nextGaussian(){return r.nextGaussian();}
	final public static double nextBiased() {  return sqrt( -1 * log(1 - nextDouble()));}//biased distribution used at walls
	final public static boolean nextBoolean(){return r.nextBoolean();}
	final public static int nextInt(int n){return r.nextInt(n);}
	final public static double abs(double x, double y, double z) {return sqrt(x * x + y * y + z * z);}	
	final public static double sq(final double x){return x*x;}
	final public static float sq(final float x){return x*x;}
	final public static double twoSig(double x){return Math.floor(x*100)/100;}
	final public static void arrayAddTo(int[] a,int[] b){for(int i=0;i<a.length;i++)a[i]+=b[i];}
	final public static void arrayAssignTo(int[] a,int x){for(int i=0;i<a.length;i++)a[i]=x;}
	final public static void arrayAssignTo(double[] d,double x){for(int i=0;i<d.length;i++)d[i]=x;}
	final public static double avg(double[] x){double res=0;for(int i=0;i<x.length;i++)res+=x[i];return res/x.length;}
	
	final public static double arrayMax(double[] a){double res=a[0];
	for(int i=0;i<a.length;i++)if(res<a[i])res=a[i];
	return res;
	}
	final public static void die(int i){System.exit(i);}
	final public static void die(int i,String st){System.err.println(st);die(i);}
	final public static double arrayMin(double[] a){double res=a[0];
	for(int i=0;i<a.length;i++)if(res>a[i])res=a[i];
	return res;
	}
	final public static double sum(final double[] x,final ArrayListInt ptrs){
		double res=0;
		for(int i=0;i<ptrs.length;i++)res+=x[ptrs.get(i)];
		return res;
	}
	final public static double avg(final double[] x,final ArrayListInt ptrs){return sum(x,ptrs)/ptrs.length;}
	final public static double variance(double[] x,ArrayListInt ptrs){
		double x2=0,ns=ptrs.length;
		for(int i=0;i<ptrs.length;i++)x2+=sq(x[ptrs.get(i)]);
		return ((x2-ns*sq(avg(x,ptrs)))/(ns-1));
	}
	final public static double variance(double[] x){
		double x2=0,ns=x.length;
		for(int i=0;i<x.length;i++)x2+=sq(x[i]);
		return ((x2-ns*sq(avg(x)))/(ns-1));
	}
	final public static void deleteFileIfExist(String fileName) {
		// A File object to represent the filename
		File f = new File(fileName);
		// Make sure the file or directory exists and isn't write protected
		if (!f.exists()) return;

		if (!f.canWrite())
			throw new IllegalArgumentException("Delete: write protected: "
					+ fileName);

		// If it is a directory, make sure it is empty
		if (f.isDirectory()) {
			String[] files = f.list();
			if (files.length > 0)
				throw new IllegalArgumentException(
						"Delete: directory not empty: " + fileName);
		}

		// Attempt to delete it
		boolean success = f.delete();

		if (!success)
			throw new IllegalArgumentException("Delete: deletion failed");
	}
	void pushArrayAt(double[] tmp, double[] w, ArrayListInt ptrs) { 
		for(int i=0;i<ptrs.length;i++){int pk=ptrs.get(i);w[pk]=tmp[pk];}
	}

	/** PDF of Maxwell boltzmann for @ a point and for a set of given parameters
	 * 
	 * @param cx
	 * @param cy
	 * @param cz
	 * @param n0
	 * @param mass
	 * @param T0
	 * @param ux
	 * @param uy
	 * @param uz
	 * @return density fmb (ie exact and absolution function multiplied by the actual number density
	 */ 
	
	final static double nfmb(final double cx,final double cy,final double cz,
			final double ux,final double uy,final double uz,
			final double n0,final double mass,final double T0
			){
		return 
		(1.2376694183279517e33*pow(mass,1.5)*n0)/
		   (exp((3.6214818480827467e22*mass*(sq(cx - ux) + sq(cy - uy) + sq(cz - uz)))/T0)
				   *pow(T0,1.5));
	}

	/** this is a generic file loader 
	 * 	it will take care of it's exceptions
	 * @param filename  name of file
	 * @return Object that was in a file
	 */
	final public static Object loadObjectFromFile(String filename){
		Object obj=null;
		FileInputStream f_in;
		try {
			f_in = new FileInputStream(filename);
			// Read object using ObjectInputStream
			ObjectInputStream obj_in = 	new ObjectInputStream (f_in);
			// Read an object
			obj = obj_in.readObject();
			if(obj==null) {throw new Exception();}
		} catch (Exception e) {
			System.err.println("readFromFile:Problem loading file "+filename);
			e.printStackTrace();
		}
		return obj;
	}
	/** my implementation of the Error Function. I copied this from:
	 * http://www.cs.princeton.edu/introcs/21function/MyMath.java.html I don't know what the copyright is 
	 * @param z
	 * @return
	 */
	
	double[] erfValues;
	static double increments = 0.01;
	
	void test() {
		
	}
	
	final public static double erf(final double z) {
		double t = 1.0 / (1.0 + 0.5 * Math.abs(z));

		// use Horner's method
		double ans = 1 - t * Math.exp( -z*z   -   1.26551223 +
				t * ( 1.00002368 +
						t * ( 0.37409196 + 
								t * ( 0.09678418 + 
										t * (-0.18628806 + 
												t * ( 0.27886807 + 
														t * (-1.13520398 + 
																t * ( 1.48851587 + 
																		t * (-0.82215223 + 
																				t * ( 0.17087277))))))))));
		if (z >= 0) return ans;
		else        return -ans;
	}
	final public static double cheapErf(final double z) {
		final double a=0.1400122886866665f;
		final double z2=(z*z);
		final double az2=a*z2;
		double ans=1.0-Math.exp(-z2*((1.273239544735162+az2)/(1+az2)));
		if (z >= 0) return sqrt(ans);
		else        return -sqrt(ans);
	}



	/** this is a generic file saver of objects. 
	 * 
	 * @param filename name of file to write to 
	 * @param obj whatever u wanna save....
	 */
	final public static void saveObject(String filename,Object obj){
//		this code is from: http://www.javacoffeebreak.com/articles/serialization/index.html
		try {
			deleteFileIfExist(filename);
			// Write to disk with FileOutputStream
			FileOutputStream f_out = new FileOutputStream(filename);
			// Write object with ObjectOutputStream
			ObjectOutputStream obj_out = new ObjectOutputStream (f_out);
			// Write object out to disk
			obj_out.writeObject(obj);
			obj_out.flush();
			System.err.print("...");//indicate we wrote things to disk. 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}


}

class State{
	double ux,uy,uz,n0,T0;
	
	public void makeRef(Parameters pars) {	
		ux=uy=uz=0;n0=pars.n0Ref;T0=pars.TRef;
	}
/** assign the VR properties
 * 
 * @param cell
 */
	public void VRmeasureAndAssign(Cell cell) {
		ux=cell.ux();uy=cell.uy();uz=cell.uz();
		n0=cell.n0();
		T0=cell.T();
	}
}
