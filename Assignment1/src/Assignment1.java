/**
 * filled in by @Casey
 */
import java.util.ArrayList;
import java.util.Collections;

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


    /**
     * Transform on one dimension
     * @param f_xy original image, but it is complexed
     * @param width
     * @param height
     * @param forward an double value, -1.0 means doing FFT, 1.0 means doing IFFT. It will be used in calculating W
     * @return F_xv (transform in one dimension)
     */
    public Complex[] fxvHelper(Complex[] f_xy, int width, int height, double forward){
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
        Complex[] F_xv_even = fxvHelper(f_xy_even, K, height, forward);
        Complex[] F_xv_odd = fxvHelper(f_xy_odd, K, height, forward);
        //conquer
        for(int x = 0; x < height; x++){
            for(int v = 0; v < K; v++){
                Complex F_xv_even_tmp = F_xv_even[x * K + v];
                Complex F_xv_odd_tmp_mul_by_w = F_xv_odd[x * K + v].mul(Complex.fromPolar(1, forward * Math.PI * v/(double)K));
                F_xv[x * width + v] = F_xv_even_tmp.plus(F_xv_odd_tmp_mul_by_w);
                F_xv[x * width + v + K] = F_xv_even_tmp.minus(F_xv_odd_tmp_mul_by_w);
            }
        }
        return F_xv;
    }

    /**
     * Transform on the other dimension
     * @param F_xv An image that has be transform into frequency domain in one dimension
     * @param width
     * @param height
     * @param forward an double value, -1.0 means doing FFT, 1.0 means doing IFFT. It will be used in calculating W
     * @return Frequency domain of the original img.
     */
    public Complex[] fuvHelper(Complex[] F_xv, int width, int height, double forward){

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
        Complex[] F_uv_even = fuvHelper(F_xv_even, width, F_xv_even.length / width, forward);
        Complex[] F_uv_odd = fuvHelper(F_xv_odd, width, F_xv_odd.length / width, forward);
        //conquer
        for(int u = 0; u < K; u++){
            for(int v = 0; v < width; v++){
                Complex F_uv_even_tmp = F_uv_even[u * width + v];
                Complex F_uv_odd_tmp_mul_by_w = F_uv_odd[u * width + v].mul(Complex.fromPolar(1, forward * Math.PI * u / (double) K));
                F_uv[u * width + v] = F_uv_even_tmp.plus(F_uv_odd_tmp_mul_by_w);
                F_uv[(u+K) * width + v] = F_uv_even_tmp.minus(F_uv_odd_tmp_mul_by_w);
            }
        }
        return F_uv;
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
        for(int x = 0; x < height; x++){
            for(int y = 0; y < width; y++){
                f_xy[x * width + y] = new Complex((img[x * width + y] & 0xff) * Math.pow(-1, x+y), 0);
            }
        }
        // prepare F_xv by first fft recursion on y axis
        Complex[] F_xv = fxvHelper(f_xy, width, height, -1.0);
        // start second recursion on x axis
        Complex[] F_uv = fuvHelper(F_xv, width, height, -1.0);
        // test if DFT and FFT are same
//        Complex[] dft = ft(img, width, height);
//        for(int i = 0; i < F_uv.length; i++){
//            if ((F_uv[i].minus(dft[i]).getNorm()<0.001)){
//                continue;
//            }
//            else {
//                System.out.println("Hehe");
//            }
//        }
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
//        for(int i = 0; i < img.length; i++){
//            if((img[i]&0xff) >0){
//                int x = i / width;
//                int y = i % width;
//                System.out.println("x:" + x + ", y:" + y + ", value:" + (img[i]&0xff));
//            }
//        }
        long endTime = System.nanoTime();
        System.out.print("run time: " + (endTime - startTime) / 1000000 + "ms");
	}

    //naive dft
    private Complex[] ft(byte[] img, int width, int height) {
        Complex[] F = new Complex[width*height];
        for (int u = 0;u<height;u++) {
            for (int v = 0;v<width;v++) {
                F[u*width + v] = new Complex();
                for (int x = 0;x<height;x++) {
                    for (int y = 0;y<width;y++) {
                        double a = (double)(img[x*width+y]&0xFF)*Math.pow(-1, x+y);
                        double b = -2.0 * Math.PI * (u*x/(double)height + v*y/(double)width);
                        Complex c = Complex.fromPolar(1, b).mul(a);
                        F[u*width+v] = F[u*width+v].plus(c);
                    }
                }
            }
        }
        return F;
    }


    //for comparison
    private Complex[] fft(byte[] img, int width, int height){
        Complex[] f_xy = new Complex[img.length];
        // transfer byte img into complex img
        for(int x = 0; x < height; x++){
            for(int y = 0; y < width; y++){
                f_xy[x * width + y] = new Complex((img[x * width + y] & 0xff) * Math.pow(-1, x+y), 0);
            }
        }
        // prepare F_xv by first fft recursion on y axis
        Complex[] F_xv = fxvHelper(f_xy, width, height, -1.0);
        // start second recursion on x axis
        return fuvHelper(F_xv, width, height, -1.0);
    }



    public void ifft(Complex[] img_ft, byte[] img, int width, int height){
        Complex[] F_uy = fxvHelper(img_ft, width, height, 1.0);
        Complex[] f_xy = fuvHelper(F_uy, width, height, 1.0);
        for(int x = 0; x < height; x++){
            for(int y = 0; y < width; y++){
                f_xy[x * width + y] = f_xy[x * width + y].div(height * width);
                f_xy[x * width + y] = f_xy[x * width + y].mul(Math.pow(-1, (x+y)));//center
                if(f_xy[x * width + y].getReal() < 0) f_xy[x * width + y].setReal(0.0);
                img[x * width + y] = (byte)(f_xy[x * width + y].getReal());
            }
        }

    }
	/**
	  * Task 2
	  * Describe the appearance of the Fourier spectrum that is obtain from the image created in this function.
	  * Why does the Fourier spectrum appear as it does? Provide a detailed explanation as a comment to this function.
	  * No need to implement any code here.
	  * 2 Marks
	  */
    /**
     * In the image created by this function, we only need to care about the x-axis because only x has variance that contributes the frequency component.
     * In the expression, the frequency of this signal is 0.125 with respect to x. Thus, if this signal is transformed into the frequency domain,
     * the image should has two points at 0.125Hz and -0.125Hz. However, the transformed image has 8 points. The following is the reason:
     *
     * There are 8 points in a period. Take the first period as an example, which x is from 0 to 7. Then the value of pixel grey level should be
     * 255, 180, 0, -180, -255, -180, 0, 180 then to 255, which is the first point of the next period. However, these values need to be
     * transformed into "byte" type value, which is an 8-bit signed type. For example, the binary code of "255" is 11111111 unsigned.
     * If it is converted into byte, the value will be changed into -1 because the complement of -1 is 11111111. The same reason to other value.
     * Thus the original value will be changed into -1, -76, 0, 76, 1, 76, 0, -76. This series will be stored into the img list.
     *
     * Next, in the fourier transform method, these value will be "and" with 0xff. Thus, the value in one period will be changed into
     * 255, 180, 0, 76, 1, 76, 0, 180. You can see the value is back to the original if they are positive. However, those negative value at
     * the very beginning such as -180 or -255 will be changed into 76 and 1 relatively. It is not a cosine wave anymore. This is the reason why we
     * will see six extra frequency components, except for the original two.
     *
     * The reason why there are 8 frequency component is that, in a period, there are 8 points in a period. It needs 8 different frequency component to
     * express it. Their frequency span is 0.125Hz.
     */
  public void changeImage(byte[] img, int width, int height) {
      for (int x = 0;x<height;x++)
        for (int y = 0;y<width;y++){
            img[x*width+y] = (byte)(255.0 * Math.cos(x / 8.0 * Math.PI));
            System.out.println(img[x * width + y] & 0xff);
        }
  }

	/**
	  * Task 3
	  * Perform second order (n=2) ButterWorth low pass filtering in the frequency domain.
	  * 5 Marks
	  */

  public void filtering(byte[] img, int width, int height, double d0) {
      int n = 2;
      Complex[] F_uv = fft(img, width, height);
      for(int u = 0; u < height; u++){
          for(int v = 0; v < width; v++){
              double D_uv = Math.sqrt(Math.pow(u - height/2, 2) + Math.pow(v - width/2, 2));
              double H_uv = 1.0 / (1.0 + Math.pow((D_uv/d0), 2 * n));
              F_uv[u * width + v] = F_uv[u * width + v].mul(H_uv);
          }
      }
      ifft(F_uv, img, width, height);
  }

	/**
	  * Task 4
	  * Apply a suitable filter to the car image to attenuate the "impulse-like" bursts in the image.
	  * 2 Mark
      * Using notch Butterworth high pass filter
	  */
  public void filtering2(byte[] img, int width, int height) {
      int d0 = 8;
      int n = 1;
      Complex[] F_uv = fft(img, width, height);
      int[] pair1 = {21, 21};
      int[] pair2 = {21, -21};
      double[] H_uv = new double[F_uv.length];
      for(int u = 0; u < height; u++){
          for(int v = 0; v < width; v++){
              double D1_uv = Math.sqrt(Math.pow(u - height/2 - pair1[0], 2) + Math.pow(v - width/2-pair1[1], 2));
              double Dm1_uv = Math.sqrt(Math.pow(u - height/2 +pair1[0], 2) + Math.pow(v - width/2+pair1[1],2));
              double D2_uv = Math.sqrt(Math.pow(u - height/2 - pair2[0], 2) + Math.pow(v - width/2-pair2[1], 2));
              double Dm2_uv = Math.sqrt(Math.pow(u - height/2 +pair2[0], 2) + Math.pow(v - width/2+pair2[1],2));
              H_uv[u * width + v] = (1.0/(1.0 + Math.pow(d0/D1_uv, 2 * n))) * (1.0/(1.0 + Math.pow(d0/Dm1_uv, 2 * n))) *
                      (1.0/(1.0 + Math.pow(d0/D2_uv, 2 * n))) * (1.0/(1.0 + Math.pow(d0/Dm2_uv, 2 * n)));
              F_uv[u * width + v] = F_uv[u * width + v].mul(H_uv[u * width + v]);
          }
      }

      ifft(F_uv, img, width, height);
  }
    //not used in this assignment
    private void medianFilter(byte[] img, int width, int height) {
        byte[] new_img = img.clone();
        int filterSize = 5;
        int halfFilterSize = filterSize / 2;
        for (int i = 0; i < new_img.length; i++) {
            int x = i / width;
            int y = i % width;
            ArrayList<Integer> effectiveAdj = new ArrayList<>();
            for (int j = -halfFilterSize; j <= halfFilterSize; j++) {
                for (int k = -halfFilterSize; k <= halfFilterSize; k++) {
                    if ((x + j) >= 0 && (x + j) < height && (y + k) >= 0 && (y + k) < width) {
                        effectiveAdj.add((new_img[(x + j) * width + (y + k)]) & 0xff);
                    }
                }
            }
            while (effectiveAdj.size() < Math.pow(filterSize, 2)) {
                effectiveAdj.add(0);
            }
            Collections.sort(effectiveAdj);
            img[x * width + y] = (byte) (effectiveAdj.get(effectiveAdj.size() / 2) & 0xff);
        }
    }
}
