import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

public class Workshop2 {
	/**
	 * Task 1: 
     * Implement the thresholding transformation. It will set 
     * all pixels with intensity below the threshold 
     * to black and other pixels to white.
     * @param img
     * @param threshold
     */
    public void thresholdTransformation(byte[] img, int threshold) {
		for(int i = 0; i < img.length; i++){
			img[i] = (img[i]& 0xff) < threshold ? (byte)0 : (byte)255;
		}
    }
   /**
     * Task 2: 
	 * Implement the negative transformation.
	 *
	 * @param img
	 */
	public void negativeTransformation(byte[] img) {
		for(int i = 0; i < img.length; i++){
			img[i] = (byte)(255 - img[i] & 0xff);
		}
	}
	/**
	 * Task 3:
     * Implement the log transformation.
     * @param img
     */
    public void logTransformation(byte[] img)
	{
		double c = 255/Math.log(256);
		//   	for(int i = 0; i < img.length; i++){
		// 	img[i] = (byte)(c * Math.log((img[i] & 0xff)+1));
		// }
		byte[] tmp = new byte[256];
		for(int i = 0; i < 256; i ++){
			tmp[i] = (byte)(c * Math.log(i+1));
		}
		for(int i = 0; i < img.length; i++){
			img[i] = tmp[(int)(img[i] & 0xff)];
		}
    }
    /**
     * Task 4:
     * Implement bit-plane slicing
     *
     * @param img
     * @param mask - between 0 to 255 in decimal
     */
	public void bitPlaneSlicing(byte img[], int mask) {
		for (int i = 0; i < img.length; i++){
			img[i] = (byte)(img[i] & mask);
		}
		//System.out.println("Using mask: " + Integer.toBinaryString(mask));
		int tmp[] = new int[256];
    }
	/**
	 * Task 5:
	 * Calculate the histogram of the image.
	 * @param img
	 * @return the histogram
	 */
	public int[] histogram (byte[] img) {
		int[] hist = new int[256];
		for (int i = 0; i < img.length; i++){
			hist[img[i]&0xff]++;
		}
		return hist;
	}

	/**
	 * Task 6:
	 * Calculate the cumulative histogram of the image.
	 * @param img
	 * @return the histogram
	 */
	public int[] cumulativeHistogram (byte[] img) {
		System.out.println("TODO: Task 6");
		int[] histogram = histogram(img);
		int[] cumHist = new int[histogram.length];

		cumHist[0] = histogram[0];
		for(int i = 1; i < histogram.length; i++){
			cumHist[i] = (histogram[i] + cumHist[i-1]);
		}
		int[] returnHist = new int[cumHist.length];
		for(int i = 0; i < cumHist.length; i ++){
			returnHist[i] = 7*cumHist[i]/img.length;
//			returnHist[i] = Math.round(7*cumHist[i]/img.length);
		}
		return returnHist;
	}

	/**
	 * Task 7:
	 * Perform histogram equalization.
	 * @param img
	 */
	public void histogramEqualization(byte img[]) {
		int[] nj = histogram(img);
		int[] tmp = new int[256];
		for(int i = 0; i < tmp.length; i++){
			for(int j = 0; j <= i; j++){
				tmp[i] += nj[j];
			}
			tmp[i] = (byte)(tmp[i] * 255 / img.length);
		}
		for(int i = 0; i < img.length; i ++){
			img[i] = (byte)tmp[(int)(img[i]&0xff)];
		}
	}

	private ArrayList<Integer> getNeighbor3by3(byte[] img, int x1, int y1, int width, int height){
		int[] adj = {x1-1, y1-1, x1-1, y1, x1-1, y1+1, x1, y1-1, x1, y1+1, x1+1, y1-1, x1+1, y1, x1+1, y1+1};//8 points stored in a 16 elements array
		ArrayList<Integer> effAdjIntensity = new ArrayList<>();
		for(int i = 0; i < adj.length; i+=2){
			int x = adj[i];
			int y = adj[i+1];
			if(height > x && x>=0 && width > y && y >= 0){//effective
				effAdjIntensity.add(img[x*width+y] & 0xff);
			}
		}
		return effAdjIntensity;
	}

	private ArrayList<Integer> getNeighborNbyN(byte[] img, int x, int y, int width, int height, int filterSize){
		int aorb = (filterSize-1)/2;
		ArrayList<Integer> effAdjIntensity = new ArrayList<Integer>();
		for(int i = -aorb; i <= aorb; i++){
			for(int j = -aorb; j <= aorb; j++){
				int x_neighbor = x+i;
				int y_neighbor = y+j;
				if(height > x_neighbor && x_neighbor >= 0 && width > y_neighbor && y_neighbor >=0){
					effAdjIntensity.add(img[x_neighbor * width + y_neighbor] & 0xff);
				}
			}
		}
		return effAdjIntensity;
	}


	private ArrayList<Integer> get4Neighbor(byte[] img, int x1, int y1, int width, int height){//if mode is true: adjacent, false: neighbor
		int[] adj = {x1-1, y1, x1, y1-1, x1, y1+1, x1+1, y1};//4 points stored in a 8 elements array  N, W, E, S
		ArrayList<Integer> effAdjIntensity = new ArrayList<>();
		for(int i = 0; i < adj.length; i+=2){
			int x = adj[i];
			int y = adj[i+1];
			if(height > x && x>=0 && width > y && y >= 0){//effective
				effAdjIntensity.add(img[x * width + y] & 0xff);
			}
		}
		return effAdjIntensity;
	}

	private ArrayList<Integer> getDigNeighbor(byte[] img, int x1, int y1, int width, int height){//if mode is true: adjacent, false: neighbor
		int[] adj = {x1-1, y1-1, x1-1, y1+1, x1+1, y1-1, x1+1, y1+1};
		ArrayList<Integer> effAdjIntensity = new ArrayList<>();
		for(int i = 0; i < adj.length; i+=2){
			int x = adj[i];
			int y = adj[i+1];
			if(height > x && x>=0 && width > y && y >= 0){//effective
				effAdjIntensity.add(img[x * width + y] & 0xff);
			}
		}
		return effAdjIntensity;
	}



	/**
	 * Homework 1 (optional):
	 * Implement the box smoothing filter.
	 *
	 * @param img the graylevel image (row major representation)
	 * @param w width of the image
	 * @param h height of the image
	 * @param filterSize the size of the filter, which is supplied by the user
	 */
	public void boxSmoothFilter(byte[] img, int w, int h, int filterSize) {
		byte[] new_img = img.clone();
		for(int i = 0; i < new_img.length; i++){
			int x = i / w;
			int y = i % w;
//			ArrayList<Integer> effAdjIntensity = getNeighbor3by3(new_img, x, y, w, h);
			ArrayList<Integer> effAdjIntensity = getNeighborNbyN(new_img, x, y, w, h, filterSize);
			int sum = 0;
			for(int eachIntensity : effAdjIntensity){
				sum += eachIntensity;
			}
			img[i] = (byte)(sum/filterSize/filterSize & 0xff);
		}
	}
	/**
	 * Homework 2 (optional):
	 * Implement the median filter. 
	 * The java function java.util.Arrays.sort(array) can be used to sort an array.
	 *
	 * @param img the graylevel image (row major representation)
	 * @param w width of the image
	 * @param h height of the image
	 * @param filterSize the size of the filter, which is supplied by the user
	 */
	public void medianFilter(byte[] img, int w, int h, int filterSize) {
		byte[] new_img = img.clone();
		for(int i = 0; i < new_img.length; i++){
			int x = i / w;
			int y = i % w;
			ArrayList<Integer> effAdjIntensity = getNeighborNbyN(new_img, x, y, w, h, filterSize);
			Collections.sort(effAdjIntensity);
			img[i] = (byte)(effAdjIntensity.get(effAdjIntensity.size()/2) & 0xff);

		}
	}
	/** 
	 * Homework 3 (optional):
	 * Implement the Laplacian filter (isotropic mask for rotations in increments of 45 deg)
	 *
	 * @param img the graylevel image (row major)
	 * @param w width of the image
	 * @param h height of the image
	 */
	public void laplacianFilter(byte img[], int w, int h) {
		byte[] new_img = img.clone();
		for(int i = 0; i < new_img.length; i++){
			int x = i / w;
			int y = i % w;
			ArrayList<Integer> effAdjIntensity = getNeighbor3by3(new_img, x, y, w, h);
			int sum = 0;
			for(int eachIntensity : effAdjIntensity){
				sum += eachIntensity;
			}
			sum = -sum + 8 * new_img[i];
			if(sum <0){
				sum = 0;
			}
			if(sum > 255){
				sum=255;
			}
			img[i] = (byte)(sum & 0xff);
		}


	}
}
