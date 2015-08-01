import java.util.ArrayList;
public class Workshop4 {

   /**
	 * Task 1:
	 * Display the gradient image.
	 *
	 * @param f the graylevel image (row major representation)
	 * @param width of the image
	 * @param height of the image
	 */
	public void gradientImage(byte[] f, int width, int height) {
		byte[] new_f = f.clone();
		for(int i = 0; i < f.length; i++){
			int x = i / width;
			int y = i % width;
			int x_minus1 = x - 1;
			int y_minus1 = y - 1;
			x_minus1 = Math.max(0, x_minus1);
			y_minus1 = Math.max(0, y_minus1);
			int partialDerivertive_x = new_f[x * width + y] - new_f[x_minus1 * width + y];
			int partialDerivertive_y = new_f[x * width + y] - new_f[x * width + y_minus1];
			double mag = Math.sqrt(partialDerivertive_x * partialDerivertive_x * 1.0 + partialDerivertive_y * partialDerivertive_y * 1.0);
			f[i] = (byte)mag;
		}
	}
	public void boxSmoothFilter(double[] img, int w, int h, int filterSize) {
		double[] new_img = img.clone();
		int halfFilterSize = filterSize / 2;
		for(int i = 0; i < new_img.length; i++){//traverse the element in image
			int x = i / w;
			int y = i % w;
			double sum = 0;
			for(int j = -halfFilterSize; j <= halfFilterSize; j++){
				for(int k = -halfFilterSize; k <= halfFilterSize; k++){
					if((x + j) >= 0 && (x + j)< h && (y + k) >= 0 && (y + k) < w){
						sum += new_img[(x + j) * w + (y + k)];
					}
				}
			}
			img[i] = sum / filterSize / filterSize;
		}
	}

	/**
	 * Task 2:
	 * Display the response map R of the harris corner detector. For simplicity, you may use a 5x5 box smoothing filter instead of Gaussian.
	 * Set all pixels to white for which R > 100000.
	 *
	 * @param f the graylevel image (row major representation)
	 * @param width of the image
	 * @param height of the image
	 */
	public void cornerResponseImage(byte[] f, int width, int height) {
		double[] fx2 = new double[f.length];
		double[] fy2 = new double[f.length];
		double[] fxy = new double[f.length];
		double[] R = new double[f.length];
		// calculate the partial derivative of respect to x and y
		for(int x = 1; x < height; x++){
			for(int y = 1; y < width; y++){
				int fx = ((int)f[x * width + y]&0xff) - ((int)f[(x-1) * width + y]&0xff);
				int fy = ((int)f[x * width + y]&0xff) - ((int)f[x * width + y - 1]&0xff);
				fx2[x * width + y] = fx * fx;
				fy2[x * width + y] = fy * fy;
				fxy[x * width + y] = fx * fy;
			}
		}

		// fx2, fy2, fxy filter by box
		boxSmoothFilter(fx2, width, height, 5);
		boxSmoothFilter(fy2, width, height, 5);
		boxSmoothFilter(fxy, width, height, 5);


		// calculate R
		for(int i = 0; i < R.length; i ++){
			R[i] = (fx2[i] * fy2[i] - fxy[i] * fxy[i]) - 0.04 * ((fx2[i] + fy2[i]) * (fx2[i] + fy2[i]));
		}

		for(int i = 0; i < R.length; i++){
			if(R[i] > 100000){
				f[i] = (byte)255;
			}
			else f[i]=(byte)0;
		}
	}

}
