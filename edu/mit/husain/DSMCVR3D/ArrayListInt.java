package edu.mit.husain.DSMCVR3D;
/**
 * (hopefully) efficient and flexible pointer arrayList of integers
 * 
 */

import java.io.*;
import static edu.mit.husain.DSMCVR3D.myUtils.*;
import static java.lang.Math.*;


final class ArrayListInt implements Serializable{
    private static final long serialVersionUID = 7L;
    transient private int[] data;
    public int length;//external size
    @SuppressWarnings("unused")
    private ArrayListInt(){}//don't want anyone to touch this constructor

    public ArrayListInt(int size) {
        data=new int[(int) ceil(size+6*sqrt(size))];
        length=0;
    }
    
    public ArrayListInt(int[] li){
        this(li.length);
    	System.arraycopy(li, 0, this.data, 0, li.length);
    	length=li.length;
    }

    public void add(int pk){ensureCapacity(length);data[length]=pk;length++;}
    
    public void reset(){if(data==null)data=new int[length];length=0;}//the check is for the case of restarting a calculation (it has to do with serialization and the transient keyword)
    public void add2(int pk){try{data[length]=pk;length++;}catch(ArrayIndexOutOfBoundsException e){data=new int[length+10];}}
    public void drop(int pk){data[pk]=data[length-1];length--;}
    public boolean dropValue(int val){
        boolean found=false;
        for(int i=0;i<length;i++){
            found=(val==get(i));
            if(found){
                drop(i);
                break;
            }
        }
        return found;
    }
    public int get(int i){return data[i];}
    

    private void ensureCapacity(int size) {
        if(size<data.length)return;

        int[] oldData=data;//save it for a while
        int newsize=10+(int)ceil(1.1*size);
        data=new int[newsize];//create a new array that's 10% bigger
        for(int i=0;i<length;i++){data[i]=oldData[i];}
        oldData=null;//to make sure we release the old data
        System.gc();//sorry but I want to make sure we are not wasteful with RAM. 
    }


    /** Choose a particle randomly from the list data and return it to caller
     * 
     * @return "pointer" to a random particle that is stored in this IntArrayList
     */
    public int chooseRandomInt() {
        return get(nextInt(length));
   
    }
    /** the standard toString Method:
     * print the data up to the limit of our data
     */
    public String toString(){
    	int[] res=trimmedArray();
    	String ans="{";
    	int i;
    	for(i=0;i<res.length-1;i++)ans+=res[i]+",";
    	if(length>0)ans+=res[res.length-1]+"}";
    	else ans+="}";
    	return ans;}
    /** this method returns a trimmed version of the array
     * @return regular int[] array trimmed to the size of data
     */
    public int[] trimmedArray(){
        int[] res=new int[length];
        for(int i=0;i<res.length;i++)res[i]=get(i);
        return res;
    }
}





