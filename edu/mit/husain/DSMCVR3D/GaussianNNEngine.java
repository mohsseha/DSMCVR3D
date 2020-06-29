package edu.mit.husain.DSMCVR3D;
import static java.lang.Math.*;
import static edu.mit.husain.DSMCVR3D.myUtils.*;


//this engline shall only be called once. 
//this class will probably be static in the sence that all cells will need one engine. 
public class GaussianNNEngine {
	final ArrayListInt[][][] mesh;
	final ArrayListInt res;//later point to the location in the cell of the particle used to quickly remove particle
	final double SD;
	final double[] x,y,z;
	final private int defaultPkPerCell;

	/** create a gassian NN Engine
	 * @param cx velocity vecotrs in x
	 * @param cy in y 
	 * @param cz in z
	 * @param pkInCell expected number of pk in the engine total
	 * @param pkPerCell how many particles per cell in velocity space. defualts to 10.  
	 */
	public GaussianNNEngine(final double[] cx,final double[] cy,final double[] cz,int pkInCell,int pkPerCellV) {
		if(pkPerCellV<1)pkPerCellV=10;
		defaultPkPerCell=pkPerCellV;
		SD = 280;//estimateSD(cx); 
		x=cx;y=cy;z=cz;
		int dim=(int)ceil(pow(pkInCell/defaultPkPerCell,.33333333));
		mesh=new ArrayListInt[dim][dim][dim];//they are obviously empty 
		res=new ArrayListInt(5*defaultPkPerCell);
	}
/** create an engine that is based on our gas positions in velocity and average pk/cell. 
 * @param gas
 * @param params
 */
	public GaussianNNEngine(Gas gas, Parameters params) {
		this(gas.Cx,gas.Cy,gas.Cz,params.pkPerCell,params.defPkVCell);
	}

//	private double estimateSD(double[] cx) {
//		int mx=min(100,cx.length);
//		double sum=0;
//		for(int i=0;i<mx;i++)sum+=cx[i]*cx[i];//assuming the mean is zero
//		return sqrt(sum/mx);
//	}

	public static GaussianNNEngine MMACreateEngine(double[] cx,double[] cy,double[] cz,int[] ptrs){
		ArrayListInt ptr=new ArrayListInt(ptrs);
		GaussianNNEngine eng=new GaussianNNEngine(cx,cy,cz,ptr.length,-1);
		eng.sort(ptr);
		return eng;
	}
	
	
	/**reset the mesh array */
	public void reset(){
		int d=dim();
		for(int i=0;i<d;i++)
			for(int j=0;j<d;j++)
				for(int k=0;k<d;k++){
					if(mesh[i][j][k]!=null)mesh[i][j][k].reset();
				}
	}

	private void sort(final ArrayListInt ptrsCell){
		reset();//before we start sorting the pks
		for(int i=0;i<ptrsCell.length;i++)add(ptrsCell.get(i));
	}
	/** add particle at location pk to the mesh so we can later retreive it's NN fast
	 * @param pk
	 */
	public void add(int pk) {
		int i=cellIndex(x[pk]),j=cellIndex(y[pk]),k=cellIndex(z[pk]);
		if(mesh[i][j][k]==null)mesh[i][j][k]=new ArrayListInt(defaultPkPerCell);
		mesh[i][j][k].add(pk);
	}
	
	/** delete old pk
	 * @param pk
	 */
	public boolean drop(int pk){
		int i=cellIndex(x[pk]),j=cellIndex(y[pk]),k=cellIndex(z[pk]);
		ArrayListInt cell=mesh[i][j][k];
		return cell.dropValue(pk);
	}

	/** find location of particle at location using the ERF fucntion
	 * @param x location
	 * @return 
	 */
	final private int cellIndex(final double x) {
		return (int)(dim()*phi(x/SD));
	}
	final private static double SQRT2=sqrt(2.0);//optimization JVM is not nearly as smart as I hoped for. 
	final private double phi(final double x){
		return 0.5*(1+cheapErf(x/SQRT2));
	}

	final public ArrayListInt NN(int offset, double eps){
	return NNpkj(ptrs.get(offset),eps);
	}
	/**return the NN iwthin a radious of eps of pkj
	 * 
	 * @param pkj
	 * @param eps
	 * @return
	 */
	final private ArrayListInt NNpkj(int pkj, double eps){
		res.reset();//now we have no results
		int i_min=cellIndex(x[pkj]-eps),j_min=cellIndex(y[pkj]-eps),k_min=cellIndex(z[pkj]-eps);
		int i_max=cellIndex(x[pkj]+eps),j_max=cellIndex(y[pkj]+eps),k_max=cellIndex(z[pkj]+eps);
		for(int i=i_min;i<=i_max;i++)
			for(int j=j_min;j<=j_max;j++)
				for(int k=k_min;k<=k_max;k++)
				{
					ArrayListInt cell=mesh[i][j][k];
					if(cell!=null){
						for(int ii=0;ii<cell.length;ii++){
							int pki=cell.get(ii);
							if(withinEps(pki,pkj,eps))res.add(pki);
						}
					}
				}

		return res;
	}
	
	/** brute force NN search
	 * @param ptrs
	 * @param pkj
	 * @param eps
	 * @return
	 */
	public ArrayListInt NN2(ArrayListInt ptrs,int pkj,double eps){
		res.reset();
		for(int i=0;i<ptrs.length;i++){
			int pki=ptrs.get(i);
			double dx=x[pki]-x[pkj],dy=y[pki]-y[pkj],dz=z[pki]-z[pkj];
			if(sqrt(dx*dx+dy*dy+dz*dz)<eps)res.add(pki);
		}
		return res;
	}
	public int[] NN2(int[] li,int pkj,double eps){
		ArrayListInt in=new ArrayListInt(li);
		return NN2(in,pkj,eps).trimmedArray();
	}
	public int[] NN2(int[] li,int[] pkli,double eps){
		int[] res=new int[10];
		for(int i:pkli)res=NN2(li,i,eps);
		return res;
	}

	private boolean withinEps(int pki, int pkj, final double eps) {
		double dx=x[pki]-x[pkj],
		dy=y[pki]-y[pkj],
		dz=z[pki]-z[pkj];
		return ((dx*dx+dy*dy+dz*dz)<(eps*eps));
	}
	public final int dim(){return mesh.length;}
	
	public double[][] pointsInMesh(int i,int j,int k){
		ArrayListInt cell=mesh[i][j][k];
		if(cell==null) {
			double[][] nully={{0.,0.,0.}};
			return nully;
		}
		double[][] points=new double[cell.length][3];
		for(int pk=0;pk<cell.length;pk++){
			int w=cell.get(pk);
			points[pk][0]=x[w];
			points[pk][1]=y[w];
			points[pk][2]=z[w];
		}
		return points;
	}
	
	private static ArrayListInt ptrs;
	public void rebuild(final ArrayListInt pkPtrs) {//trivial wrapper for debuggin reasons 
		ptrs=pkPtrs;
		sort(pkPtrs);
	}
	public void cache(int offset) {
		drop(ptrs.get(offset));//arugment was in offsets 
	}
	public void update(int offA) {
		add(ptrs.get(offA));//again offsets -> address wrapper
	}
	
}
class SlowPoint{

    final double[] x;
    SlowPoint(double x1,double y,double z){x=new double[3];
    x[0]=x1;
    x[1]=y;
    x[2]=z;
    }
    SlowPoint(SlowPoint pt){this(pt.x[0],pt.x[1],pt.x[2]);}
    void assign(SlowPoint pt){x[0]=pt.x[0];x[1]=pt.x[1];x[2]=pt.x[2];}
    public void assignVs(double xx,double yy,double zz){x[0]=xx;x[1]=yy;x[2]=zz;}
    public String toString(){return "{"+x[0]+","+x[1]+","+x[2]+"}";}
    public void assignVsAtPos(double[] cx, double[] cy, double[] cz, int i) {
            cx[i]=x[0];cy[i]=x[1];cz[i]=x[2];
    }

    public double[] toArray(){return x;}//used for debuging

}

