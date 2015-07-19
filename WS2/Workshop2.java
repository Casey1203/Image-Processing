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
    		/*if((img[i]& 0xff)<threshold){
    			img[i] = (byte)0;
    		}
    		else{
    			img[i] = (byte)255;
    		}*/
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
    public void logTransformation(byte[] img) {
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
     * Implement the bit-plane slicing
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
        
     	System.out.println("TODO: implement Task 4");
    }
	/**
	 * Task 5:
	 * Calculate the histogram of the image.
	 * @param img
	 * @return the histogram
	 */
	public int[] histogram (byte[] img) {
		System.out.println("TODO: implement Task 5");
		int[] hist = new int[256];
		for (int i = 0; i < img.length; i++){
			hist[img[i]&0xff]++;
		}
		return hist;
	}
	/**
	 * Task 6:
	 * Perform histogram equalization.
	 * @param img
	 */
	public void histogramEqualization(byte img[]) {
		System.out.println("TODO: implement Task 6");
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
	/**
	 * Homework 1:
	 * Implement the box smoothing filter.
	 *
	 * @param img the graylevel image (row major representation)
	 * @param w width of the image
	 * @param h height of the image
	 * @param filterSize the size of the filter, which is supplied by the user
	 */
	public void boxSmoothFilter(byte[] img, int w, int h, int filterSize) {
		System.out.println("TODO: implement Homework 1");
	}
	/**
	 * Homework 2:
	 * Implement the median filter. 
	 * The java function java.util.Arrays.sort(array) can be used to sort an array.
	 *
	 * @param img the graylevel image (row major representation)
	 * @param w width of the image
	 * @param h height of the image
	 * @param filterSize the size of the filter, which is supplied by the user
	 */
	public void medianFilter(byte[] img, int w, int h, int filterSize) {
		System.out.println("TODO: implement Homework 2");
	}
	/** 
	 * Homework 3:
	 * Implement the Laplacian filter (isotropic mask for rotations in increments of 45 deg)
	 *
	 * @param img the graylevel image (row major)
	 * @param w width of the image
	 * @param h height of the image
	 */
	public void laplacianFilter(byte img[], int w, int h) {
		System.out.println("TODO: implement Homework 3");
	}
}
