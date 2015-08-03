/**
  * Assignment 1
  * You may only edit and submit this file to moodle before the deadline. 
  */
public class Assignment1 {

    private void logTransformation(byte[] img) {
        double c = 255 / Math.log(256);
        byte[] tmp = new byte[256];
        for (int i=0;i<256;i++)
            tmp[i] = (byte)(c*Math.log(i+1));
        for (int i=0;i<img.length;i++)
            img[i] = tmp[(int)(img[i] & 0xFF)];
    }


    public Complex[] ftHelper(Complex[] F_xv, int width, int height){

        if(height <= 1){
            return F_xv;
        }
        int K = height / 2;


        Complex[] F_xv_even = new Complex[K * width];
        Complex[] F_xv_odd = new Complex[K * width];
        Complex[] F_uv = new Complex[F_xv.length];

        // calculate f(x)_even and f(x)_odd
        for(int x = 0; x < K; x++){
            for(int v = 0; v < width; v++){
                F_xv_even[x * width + v] = F_xv[2 * x * width + v];
                F_xv_odd[x * width + v] = F_xv[(2 * x + 1) * width + v];
            }
        }
        //divide
        Complex[] F_uv_even = ftHelper(F_xv_even, width, F_xv_even.length / width);
        Complex[] F_uv_odd = ftHelper(F_xv_odd, width, F_xv_odd.length / width);
        //conquer
        for(int u = 0; u < K; u++){
            for(int v = 0; v < width; v++){
                Complex F_uv_even_tmp = F_uv_even[u * width + v];
                Complex F_uv_odd_tmp_mul_by_w = F_uv_odd[u * width + v].mul(Complex.fromPolar(1, -1.0 * Math.PI * u / (double) K));
                F_uv[u * width + v] = F_uv_even_tmp.plus(F_uv_odd_tmp_mul_by_w);
                F_uv[(u+K) * width + v] = F_uv_even_tmp.minus(F_uv_odd_tmp_mul_by_w);
            }
        }
        return F_uv;
    }

    public Complex[] ftxvHelper(Complex[] f_xy, int width, int height){
        if(width <= 1){
            return f_xy;
        }
        int K = width / 2;
        // with respect to y
        Complex[] f_xy_even = new Complex[K * height];
        Complex[] f_xy_odd = new Complex[K * height];
        Complex[] F_xv = new Complex[f_xy.length];
        // calculate f(y)_even and f(y)_odd
        for(int y = 0; y < K; y++){
            for(int x = 0; x < height; x++){
                f_xy_even[x * K + y] = f_xy[x * width + 2 * y];
                f_xy_odd[x * K + y] = f_xy[x * width + 2 * y + 1];
            }
        }
        //divide
        Complex[] F_xv_even = ftxvHelper(f_xy_even, K, height);
        Complex[] F_xv_odd = ftxvHelper(f_xy_odd, K, height);
        //conquer
        for(int x = 0; x < height; x++){
            for(int v = 0; v < K; v++){
                Complex F_xv_even_tmp = F_xv_even[x * K + v];
                Complex F_xv_odd_tmp_mul_by_w = F_xv_odd[x * K + v].mul(Complex.fromPolar(1, -1.0 * Math.PI * v/(double)K));
                F_xv[x * width + v] = F_xv_even_tmp.plus(F_xv_odd_tmp_mul_by_w);
                F_xv[x * width + v + K] = F_xv_even_tmp.minus(F_xv_odd_tmp_mul_by_w);
            }
        }
        return F_xv;
    }

    /**
	* Task 1
	* Implement the Fast Fourier Transform (FFT) and display the Fourier spectrum.
	* Results should be equivalent to those obtained in the Workshop 3.
	* The implementation details of the FFT can be obtained from section 4.11 of chapter 4.
	* 6 Marks
	*/

	public void fourierTransform(byte[] img, int width, int height) {
        long startTime = System.nanoTime();
        Complex[] f_xy = new Complex[img.length];
        // transfer byte img into complex img
        for(int i = 0; i < img.length; i++){
            f_xy[i] = new Complex((int)img[i] & 0xff, 0);
        }
        // prepare F_xv by first fft recursion on y axis
        Complex[] F_xv = ftxvHelper(f_xy, width, height);
        // start second recursion on x axis
        Complex[] F_uv = ftHelper(F_xv, width, height);

        // scale processing
        double max = Double.NEGATIVE_INFINITY;
        for(int i = 0; i < F_uv.length; i++) {
            max = Math.max(F_uv[i].getNorm(), max);
        }
        double c = 255.0 / max;
        for(int i = 0; i < img.length; i++){
            img[i] = (byte)(c * F_uv[i].getNorm());
        }
        logTransformation(img);

        long endTime = System.nanoTime();
        System.out.print("run time: " + (endTime - startTime) / 1000000 + "ms");
	}

	/**
	  * Task 2
	  * Describe the appearance of the Fourier spectrum that is obtain from the image created in this function.
	  * Why does the Fourier spectrum appear as it does? Provide a detailed explanation as a comment to this function.
	  * No need to implement any code here.
	  * 2 Marks
	  */
  public void changeImage(byte[] img, int width, int height) {
	for (int x = 0;x<height;x++)
	  for (int y = 0;y<width;y++)
		img[x*width+y] = (byte)(255.0 * Math.cos(x / 4.0 * Math.PI));
  }

	/**
	  * Task 3
	  * Perform second order (n=2) ButterWorth low pass filtering in the frequency domain.
	  * 5 Marks
	  */
  public void filtering(byte[] img, int width, int height, double d0) {
  }

	/**
	  * Task 4
	  * Apply a suitable filter to the car image to attenuate the "impulse-like" bursts in the image.
	  * 2 Mark
	  */
  public void filtering2(byte[] img, int width, int height) {
  }

}
