import java.lang.invoke.ConstantCallSite;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Queue;


/*
 * Compile: javac *.java
 * Run: java Workshop1UI
 */
public class Workshop1 {
	/**
	 * Task 1:
	 * Retrieve the intensity value at location (x, y) of the image and return it
	 * Note: 
	 * - the 2D image is stored as an 8bit, 1D, row-major array of type byte
	 * - Note that byte is a signed data type in Java
	 * @param img in row major format
	 * @param x coordinate
	 * @param y coordinate
	 * @param width of img
	 * @param height of img
	 * @return the intensity value at (x, y) if within bounds, -1 otherwise
	 */
	final int V = 127;
	
	class AdvancedIntArray{
		int[] pos;
		//ArrayList<int[]> path;
		AdvancedIntArray from;
		
		public AdvancedIntArray(){
			this.pos = null;
			//this.path = null;
			this.from = null;
		}
		
		public AdvancedIntArray(int[] pos, AdvancedIntArray from){
			this.pos = pos;
			this.from = from;
		}
		
		public AdvancedIntArray getFrom(){
			return this.from;
		}
		
		public int[] getPos(){
			return this.pos;
		}
		
		public boolean isFromEmpty(){
			if (this.from == null){
				return true;
			}
			else {
				return false;
			}
		}
	}
	
	
	public int getIntensityValue(byte[] img, int x, int y, int width, int height) {
		System.out.println("TODO: implement Task 1");
		return -1;
	}
	/**
	 * Task 2:
	 * Retrieve the intensity value that occurs most often in the image
	 * @param img 
	 * @return the intensity value that occurs most often in the image
	 */
	public int getMostFrequentIntensityValue(byte[] img) {
		System.out.println("TODO: implement Task 2");
		return -1;
	}
	/**
	 * Task 3:
	 * Sets the pixels N8(row, column) to white.
	 * @param img with four neigbors of (row, column) set to white
	 * @param x coordinate
	 * @param y coordinate
	 * @param width of img
	 * @param height of img
	 */


	public void setEightNeighborsToWhite(byte[] img, int x, int y, int width, int height) {
		System.out.println("TODO: implement Task 3");
		
		int curLocation = x * width + y;
		img[curLocation-1] = (byte)255;
		img[curLocation+1] = (byte)255;
		img[curLocation-width] = (byte)255;
		img[curLocation+width] = (byte)255;
		img[curLocation-width-1] = (byte)255;
		img[curLocation-width+1] = (byte)255;
		img[curLocation+width-1] = (byte)255;
		img[curLocation+width+1] = (byte)255;

	}
	/**
	 * Task 4:
	 * Calculates the d4 distance between (x1, y1) and (x2, y2)
	 * @param img that will be unchanged
	 * @param x1 
	 * @param y1
	 * @param x2
	 * @param y2
	 * @param width of img
	 * @param height of img
	 * @return the d4 distance between (x1, y1) and (x2, y2) 
	 */
	public int getD4Distance(byte[] img, int x1, int y1, int x2, int y2, int width, int height) {
		System.out.println("TODO: implement Task 4");
		return Math.abs(x1-x2) + Math.abs(y1-y2);
	}
	/**
	 * Homework 1: 
	 * Marks the shortest m-path with white intensity values. Let V = {0, ..., 127}.
	 * This task was developed to challenge yourself. Can you find a shortest m-path quickly?
	 * @param img with the shortest m-path set to white with V = {0, ..., 127}.
	 * @param x1 
	 * @param y1
	 * @param x2
	 * @param y2
	 * @param width of img
	 * @param height of img
	 */
	private ArrayList<int[]> findEightAdjacent(byte[] img, int x1, int y1, int width, int height, boolean mode){//if mode is true: adjacent, false: neighbor
		int[] adj = {x1-1, y1-1, x1-1, y1, x1-1, y1+1, x1, y1-1, x1, y1+1, x1+1, y1-1, x1+1, y1, x1+1, y1+1};//8 points stored in a 16 elements array
		ArrayList<int[]> effAdj = new ArrayList<int[]>();
		for(int i = 0; i < adj.length; i+=2){
			int x = adj[i];
			int y = adj[i+1];
			if(height > x && x>=0 && width > y && y >= 0){//effective
				int[] pos = {x,y};
				if(mode == false){//just neighbor
					effAdj.add(pos);
				}
				else if((img[x * width + y] & 0xff) <= V){
					effAdj.add(pos);
				}
			}
		}
		return effAdj;
	}
	
	private ArrayList<int[]> findFourAdjacent(byte[] img, int x1, int y1, int width, int height, boolean mode){//if mode is true: adjacent, false: neighbor
		int[] adj = {x1-1, y1, x1, y1-1, x1, y1+1, x1+1, y1};//4 points stored in a 8 elements array  N, W, E, S
		ArrayList<int[]> effAdj = new ArrayList<int[]>();
		for(int i = 0; i < adj.length; i+=2){
			int x = adj[i];
			int y = adj[i+1];
			if(height > x && x>=0 && width > y && y >= 0){//effective
				int[] pos = {x,y};
				if(mode == false){//just neighbor
					effAdj.add(pos);
				}
				else if((img[x * width + y] & 0xff) <= V){
					effAdj.add(pos);
				}
			}
		}
		return effAdj;
	}
	private ArrayList<int[]> findDNeighbor(byte[] img, int x1, int y1, int width, int height){
		int[] adj = {x1-1, y1-1, x1-1, y1+1, x1+1, y1-1, x1+1, y1+1};//4 points stored in a 8 elements array 
		ArrayList<int[]> effAdj = new ArrayList<int[]>();
		for(int i = 0; i < adj.length; i+=2){
			int x = adj[i];
			int y = adj[i+1];
			if(height > x && x>=0 && width > y && y >= 0){//effective
				int[] pos = {x,y};
				effAdj.add(pos);
			}
		}
		return effAdj;
	}
	private boolean isInDiag(int x1, int y1, int[] pos) {
		if (pos[0] != x1){
			if (pos[1] != y1){
				return true;
			}
			else {
				return false;
			}
		}
		else
			return false;
	}
	
	private ArrayList<int[]> findMadadjacent(byte[] img, int x1, int y1, int width, int height){
		ArrayList<int[]> eightAdjacent = findEightAdjacent(img, x1, y1, width, height, false);
		ArrayList<int[]> fourAdjacent = findFourAdjacent(img, x1, y1, width, height, false);//N4(p)
		
		//ArrayList<int[]> dNeighbor = findDNeighbor(img, x1, y1, width, height);//Nd(p)
		ArrayList<int[]> mSet = new ArrayList<int[]>();
		for (int i =0; i < eightAdjacent.size(); i++){
			int[] pos = eightAdjacent.get(i);//q
			if((img[pos[0] * width + pos[1]] & 0xff) > V)
				continue;
			if(pos[0] == x1 || pos[1] == y1){
				mSet.add(pos);
			}
			else if(isInDiag(x1, y1, pos)){
				if(pos[0] < x1){
					if(pos[1] < y1){
						if(img[(x1-1) * width + y1] >V && img[x1 * width + y1-1]>V){
							mSet.add(pos);
						}
					}
					else if(pos[1] > y1){
						if(img[(x1-1) * width + y1]> V && img[x1 * width + y1+1]>V){
							mSet.add(pos);
						}
					}
				}
				else if(pos[0] > x1){
					if(pos[1] < y1){
						if(img[(x1+1) * width + y1] >V && img[x1 * width + y1 - 1]>V){
							mSet.add(pos);
						}
					}
					else if(pos[1] > y1){
						if(img[(x1+1) * width + y1]> V && img[x1 * width + y1 - 1]>V){
							mSet.add(pos);
						}
					}
				}
				
			}
		}
		return mSet;
	}
	private boolean isContain(boolean[] set, int[] pos, int width){
		return set[pos[0] * width + pos[1]];
	}
	private boolean goalTest(int[] curPos, int[] goal){
		if(Arrays.equals(curPos, goal)){
			return true;
		}
		else
			return false;
	}
	
	private void setWhite(byte[] img, AdvancedIntArray currentCondition, int width, int height) {
		if(!currentCondition.isFromEmpty()){
			img[currentCondition.getPos()[0] * width + currentCondition.getPos()[1]] = (byte)0xFF;
			System.out.println(Arrays.toString(currentCondition.getPos()));
			setWhite(img, currentCondition.getFrom(), width, height);
		}
	}
	
	public void setShortestMPathToWhite(byte[] img, int x1, int y1, int x2, int y2, int width, int height) {
		// System.out.println("TODO: implement Homework 1");
		if((img[x1 * width + y1] & 0xff) > V){
			System.out.println("no m-Path from "+ x1 + "," + y1);
			return;
		}
		System.out.println("x1, y1 = "+ x1 + "," + y1);
		System.out.println("x2, y2 = "+ x2 + "," + y2);
		int[] goal = {x2, y2};
		Queue<AdvancedIntArray> frontier = new LinkedList<Workshop1.AdvancedIntArray>();
	 	//ArrayList<int[]> exploredSet = new ArrayList<int[]>();
	 	boolean[] exploredSet = new boolean[width * height];
	 	for(int i = 0; i < exploredSet.length; i++){
	 		exploredSet[i] = false;//init
	 	}
	 	int[] cur = {x1, y1};
	 	//exploredSet.add(cur);
	 	exploredSet[x1 * width + y1] = true;
	 	AdvancedIntArray ai = new AdvancedIntArray(cur, new AdvancedIntArray());
	 	frontier.offer(ai);
	 	int count = 0;
	 	while(!frontier.isEmpty()){
	 		count++;
	 		AdvancedIntArray currentCondition = frontier.poll();


	 		//System.out.println(Arrays.toString(currentCondition.getPos()));
	 		//System.out.println("test " + count + " point");
	 		if(goalTest(currentCondition.getPos(), goal)){
	 			setWhite(img, currentCondition, width, height);
	 			System.out.println("after" + count + " find the goal");
	 			return;
	 		}
	 		else {
	 			ArrayList<int[]> mAdjList = findMadadjacent(img, currentCondition.getPos()[0], currentCondition.getPos()[1], width, height);
	 			for(int i = 0; i < mAdjList.size(); i++){
	 				if(!isContain(exploredSet, mAdjList.get(i), width)){
	 					//exploredSet.add(mAdjList.get(i));
	 					exploredSet[mAdjList.get(i)[0] * width + mAdjList.get(i)[1]] = true;
	 					frontier.offer(new AdvancedIntArray(mAdjList.get(i), currentCondition));
	 				}
	 			}
	 		}
	 	}
	}
	
	
	//  public static void main(String[] args){
	//  	Workshop1 ws = new Workshop1();
	//  	byte[] img = {3,1,2,1,2,2,0,2,1,2,1,1,1,0,1,2};
	//  	int x = 3, y = 0;
	//  	//System.out.println(Arrays.toString(ws.findMadadjacent(img, x, y, 4, 4).get(1)));
	//  	ws.setShortestMPathToWhite(img, x, y, 0, 3, 4, 4);	
	// }
	
	
	
	
	
}
