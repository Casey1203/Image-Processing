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




    public Complex[] fxvHelper(Complex[] f_xy, int width, int height){
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
        Complex[] F_xv_even = fxvHelper(f_xy_even, K, height);
        Complex[] F_xv_odd = fxvHelper(f_xy_odd, K, height);
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
    public Complex[] fuyHelper(Complex[] F_uv, int width, int height){
        if(width <= 1){
            return F_uv;
        }
        int K = width / 2;
        Complex[] F_uv_even = new Complex[K * height];
        Complex[] F_uv_odd = new Complex[K * height];
        Complex[] F_uy = new Complex[F_uv.length];
        for(int u = 0; u < height; u++){
            for(int v = 0; v < K; v++){
                F_uv_even[u * K + v] = F_uv[u * width + 2 * v];
                F_uv_odd[u * K + v] = F_uv[u * width + 2 * v + 1];
            }
        }
        //divide
        Complex[] F_uy_even = fuyHelper(F_uv_even, K, height);
        Complex[] F_uy_odd = fuyHelper(F_uv_odd, K, height);
        //conquer
        for(int u = 0; u < height; u++){
            for(int y = 0; y < K; y++){
                Complex F_uy_even_tmp = F_uy_even[u * K + y];
                Complex F_uy_odd_tmp_my_by_w = F_uy_odd[u * K + y].mul(Complex.fromPolar(1, 1.0 * Math.PI * y/(double)K));
                F_uy[u * width + y] = F_uy_even_tmp.plus(F_uy_odd_tmp_my_by_w);
                F_uy[u * width + y + K] = F_uy_even_tmp.minus(F_uy_odd_tmp_my_by_w);
            }
        }
        return F_uy;
    }

    public Complex[] fuvHelper(Complex[] F_xv, int width, int height){

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
        Complex[] F_uv_even = fuvHelper(F_xv_even, width, F_xv_even.length / width);
        Complex[] F_uv_odd = fuvHelper(F_xv_odd, width, F_xv_odd.length / width);
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

    public Complex[] fxyHelper(Complex[] F_uy, int width, int height){
        if(height <= 1){
            return F_uy;
        }
        int K = height / 2;
        Complex[] F_uy_even = new Complex[K * width];
        Complex[] F_uy_odd = new Complex[K * width];
        Complex[] f_xy = new Complex[F_uy.length];
        // calculate F(u)_even and F(u)_odd
        for(int u = 0; u < K; u++){
            for(int y = 0; y < width; y++){
                F_uy_even[u * width + y] = F_uy[2 * u * width + y];
                F_uy_odd[u * width + y] = F_uy[(2 * u + 1) * width + y];
            }
        }
        //divide
        Complex[] f_xy_even = fxyHelper(F_uy_even, width, K);
        Complex[] f_xy_odd = fxyHelper(F_uy_odd, width, K);
        //conquer
        for(int x = 0; x < K; x++){
            for(int y = 0; y < width; y++){
                Complex f_xy_even_tmp = f_xy_even[x * width + y];
                Complex f_xy_odd_tmp_mul_by_w = f_xy_odd[x * width + y].mul(Complex.fromPolar(1, 1.0 * Math.PI * x/(double)K));
                f_xy[x * width + y] = f_xy_even_tmp.plus(f_xy_odd_tmp_mul_by_w);
                f_xy[(x+K) * width + y] = f_xy_even_tmp.minus(f_xy_odd_tmp_mul_by_w);
            }
        }
        return f_xy;
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
            f_xy[i] = new Complex((img[i] & 0xff) * Math.pow(-1, i), 0);
        }
        // prepare F_xv by first fft recursion on y axis

        Complex[] F_xv = fxvHelper(f_xy, width, height);
        // start second recursion on x axis
        for(int i = 0; i < F_xv.length; i++){
            F_xv[i] = F_xv[i].mul(Math.pow(-1, i / width));
        }
        Complex[] F_uv = fuvHelper(F_xv, width, height);

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

    private Complex[] fft(byte[] img, int width, int height){
        Complex[] f_xy = new Complex[img.length];
        // transfer byte img into complex img
        for(int i = 0; i < img.length; i++){
            f_xy[i] = new Complex((img[i] & 0xff) * Math.pow(-1, i), 0);
        }
        // prepare F_xv by first fft recursion on y axis
        Complex[] F_xv = fxvHelper(f_xy, width, height);
        // start second recursion on x axis
//        for(int i = 0; i < F_xv.length; i++){
//            F_xv[i] = F_xv[i].mul(Math.pow(-1, i / width));
//        }
        return fuvHelper(F_xv, width, height);
    }



    public void ifft(Complex[] img_ft, byte[] img, int width, int height){

        Complex[] F_uy = fuyHelper(img_ft, width, height);
        Complex[] f_xy = fxyHelper(F_uy, width, height);
        for(int i = 0; i < f_xy.length; i++){
            f_xy[i] = f_xy[i].div(f_xy.length);
            if(f_xy[i].getReal() < 0) f_xy[i].setReal(0.0);
            img[i] = (byte)(f_xy[i].getReal());
        }
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
      Complex[] F_uv = fft(img, width, height);
      for(int u = 0; u < height/2; u++){
          for(int v = 0; v < width/2; v++){
              double D_uv = Math.sqrt(u * u + v * v);
              double H_uv = 1.0 / (1.0 + (D_uv/d0) * (D_uv/d0) * (D_uv/d0) * (D_uv/d0));
              F_uv[u * width + v] = F_uv[u * width + v].mul(H_uv);
          }
      }
      for(int u = height-1; u >= height/2; u--){
          for(int v = width-1; v >= width/2; v--){
              double D_uv = Math.sqrt(Math.pow(u-height,2) + Math.pow((v-width),2));
              double H_uv = 1.0 / (1.0 + (D_uv/d0) * (D_uv/d0) * (D_uv/d0) * (D_uv/d0));
              F_uv[u * width + v] = F_uv[u * width + v].mul(H_uv);
          }
      }
      for(int u = 0; u < height/2; u++){
          for(int v = width-1; v>=width/2; v--){
              double D_uv = Math.sqrt(Math.pow(u, 2) + Math.pow((v-width), 2));
              double H_uv = 1.0 / (1.0 + (D_uv/d0) * (D_uv/d0) * (D_uv/d0) * (D_uv/d0));
              F_uv[u * width + v] = F_uv[u * width + v].mul(H_uv);
          }
      }
      for(int u = height-1; u >= height/2; u--){
          for(int v = 0; v<width/2; v++){
              double D_uv = Math.sqrt(Math.pow(u-height, 2) + Math.pow(v, 2));
              double H_uv = 1.0 / (1.0 + (D_uv/d0) * (D_uv/d0) * (D_uv/d0) * (D_uv/d0));
              F_uv[u * width + v] = F_uv[u * width + v].mul(H_uv);
          }
      }
      ifft(F_uv, img, width, height);
  }

	/**
	  * Task 4
	  * Apply a suitable filter to the car image to attenuate the "impulse-like" bursts in the image.
	  * 2 Mark
	  */
  public void filtering2(byte[] img, int width, int height) {
  }

}
