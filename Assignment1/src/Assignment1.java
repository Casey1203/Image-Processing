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
	* Task 1
	* Implement the Fast Fourier Transform (FFT) and display the Fourier spectrum.
	* Results should be equivalent to those obtained in the Workshop 3.
	* The implementation details of the FFT can be obtained from section 4.11 of chapter 4.
	* 6 Marks
	*/

	public void fourierTransform(byte[] img, int width, int height) {
        byte[] new_img = img.clone();
        Complex[] F_xv = new Complex[img.length];
        Complex[] F_uv = new Complex[img.length];
        int M = width;
        int K = M / 2;
//        for(int v = 0; v < width; v++){
//            for(int x = 0; x < height; x++){
//                Complex sum = new Complex();
//                for(int y = 0; y < width; y++){
//                    double mag = new_img[x * width + y];
//                    double angle = -2 * Math.PI * (v * y / (double)width);
//                    Complex tmp = Complex.fromPolar(mag, angle);
//                    sum = sum.plus(tmp);
//                }
//                F_xv[x * width + v] = sum;
//            }
//        }
        for(int x = 0; x < height; x++){
            for(int v = 0; v < width; v++){
                Complex sum = new Complex();
                for(int y = 0; y < width; y++){
                    double mag = new_img[x * width + y] * Math.pow(-1, y);
                    double angle = -2 * Math.PI * v * y /(double)width;
                    Complex tmp = Complex.fromPolar(mag, angle);
                    sum = sum.plus(tmp);
                }
                F_xv[x * width + v] = sum;
            }
        }
//        for(int u = 0; u < height; u++){
//            for(int v = 0; v < width; v++){
//                Complex sum = new Complex();
//                for(int x = 0; x < height; x++){
//                    Complex mag = F_xv[x * width + v];
//                    double angle = -2 * Math.PI * u * x / (double)height;
//                    Complex eular = Complex.fromPolar(1 * Math.pow(-1, x), angle);
//                    sum = sum.plus(mag.mul(eular));
//                }
//                F_uv[u * width + v] = sum;
//            }
//        }
        for(int u = 0; u < K; u++){
            for(int v = 0; v < height; v++){
                F_uv[u * width + v] = new Complex();
                F_uv[(u+K) * width + v] = new Complex();

                for(int x = 0; x < K ; x++){
                    // even part
                    Complex mag_even = F_xv[2 * x * width + v];
                    double angle_even = -2 * Math.PI * (u * x / (double)K);
                    Complex tmp_even = Complex.fromPolar(1, angle_even).mul(mag_even);
                    F_uv[u * width + v] = F_uv[u * width + v].plus(tmp_even);
                    F_uv[(u+K) * width + v] = F_uv[(u+K) * width + v].plus(tmp_even);

                    // odd part
                    Complex mag_odd = F_xv[(2 * x + 1) * width + v];
                    double angle_odd = -2 * Math.PI * ((u * x + 0.5 * u) / (double)K);
                    Complex tmp_odd = Complex.fromPolar(-1, angle_odd).mul(mag_odd);
                    F_uv[u * width + v] = F_uv[u * width + v].plus(tmp_odd);
                    F_uv[(u+K) * width + v] = F_uv[(u+K) * width + v].minus(tmp_odd);

                }
            }
        }
        double max = Double.NEGATIVE_INFINITY;
        for(int i = 0; i < F_uv.length; i++) {
            max = Math.max(F_uv[i].getNorm(), max);
        }
        double c = 255.0 / max;
        for(int i = 0; i < img.length; i++){
            img[i] = (byte)(c * F_uv[i].getNorm());
        }
//        logTransformation(img);

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
