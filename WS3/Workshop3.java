import javax.activation.MailcapCommandMap;
import java.util.Collection;

public class Workshop3 {
	/**
	 * Log transformation. You may make use of this function when solving Task 1.
	 *
	 * @param img the graylevel image (row major representation)
	 */
    private void logTransformation(byte[] img) {
    	double c = 255 / Math.log(256);
		byte[] tmp = new byte[256];
		for (int i=0;i<256;i++)
			tmp[i] = (byte)(c*Math.log(i+1));
		for (int i=0;i<img.length;i++)
			img[i] = tmp[(int)(img[i] & 0xFF)];
	}

	public Complex[] ft(byte[] img, int width, int height){
		Complex[] ft = new Complex[img.length];
		for(int u = 0; u < height; u++){
			for(int v = 0; v < width; v++){
				Complex sum = new Complex();
				for(int x = 0; x < height; x++){
					for(int y = 0; y < width; y++){
						double mag = (img[x * width + y] & 0xff);
						double angle = -2 * Math.PI * (u * x/(double)height + v * y/(double)width);
						Complex tmp = Complex.fromPolar(mag, angle);
						sum = sum.plus(tmp);
					}
				}
				ft[u * width + v] = sum;
			}
		}
		return ft;
	}


	/**
	 * Task 1:
	 * Implement the Fourier transform. Make sure all intensity values are in the range 0 .. 255.
	 *
	 * @param img the graylevel image (row major representation)
	 * @param width of the image
	 * @param height of the image
	 */
    public void fourierTransform(byte[] img, int width, int height) {
//		byte[] new_img = img.clone();
//		for(int u = 0; u < height; u++){
//			for(int v = 0; v < width; v++){
//				double sum = 0;
//				for(int x = 0; x < height; x++){
//					for(int y = 0; y < width; y++){
//						double angle = 2*Math.PI*(u*x/(double)height+v*y/(double)width);
//						double real = Math.cos(angle);
//						double imag = -Math.sin(angle);
//						sum += (double)(new_img[x * width + y] & 0xff) * Math.sqrt(real * real + imag * imag)/Math.sqrt(width * height);
//
//					}
//				}
//
//				img[u * width + v] = (byte)sum;
//			}
//		}

		Complex[] ft = new Complex[img.length];
		for(int u = 0; u < height; u++){
			for(int v = 0; v < width; v++){
				Complex sum = new Complex();
				for(int x = 0; x < height; x++){
					for(int y = 0; y < width; y++){
						double mag = (img[x * width + y] & 0xff) * Math.pow(-1, x+y);
						double angle = -2 * Math.PI * (u * x/(double)height + v * y/(double)width);
						Complex tmp = Complex.fromPolar(mag, angle);
						sum = sum.plus(tmp);
					}
				}
				ft[u * width + v] = sum;
				double norm = sum.getNorm()/Math.sqrt(width * height);
				norm = Math.max(0, norm);
				norm = Math.min(norm, 255);
//				System.out.println(norm);
				img[u * width + v] = (byte)norm;
			}
		}
//		logTransformation(img);
    }

	public void idft(Complex[] img_ft, byte[] img, int width, int height){
		for(int x = 0; x < height; x++){
			for(int y = 0; y < width; y++){
				Complex sum = new Complex();
				for(int u = 0; u < height; u++){
					for(int v = 0; v < width; v++){
						Complex tmp = new Complex(Math.cos(2 * Math.PI * (u * x / (double) height + (v * y) / (double) width)), Math.sin(2 * Math.PI * (u * x / (double) height + (v * y) / (double) width)));
						sum = sum.plus(img_ft[u * width + v].mul(tmp));
					}
				}

				double norm = sum.getNorm()/(width * height);
//				System.out.println(norm);
				norm = Math.max(0, norm);
				norm = Math.min(norm, 255);

				img[x * width + y] = (byte)norm;
			}
		}
	}

   /**
	 * Task 2:
	 * Change the image to test your Fourier transform implementation. 
	 *
	 * @param img the graylevel image (row major representation)
	 * @param width of the image
	 * @param height of the image
	 */
    public void changeImage(byte[] img, int width, int height) {
		for(int i = 0; i < height; i++){
			for(int j = 0; j < width; j++){
				if(j > width/2-20 && j <width/2+20 && i > height/2-10 && i < height/2+10){
					img[i * width + j] = (byte)255;
				}
				else{
					img[i * width + j] = (byte)0;
				}
			}
		}
	}


   /**
	 * Task 3:
	 * Implement a simple filter, e.g. one that sets the center of the transform to 0. 
	 *
	 * @param img the graylevel image (row major representation)
	 * @param width of the image
	 * @param height of the image
	 */
	public void filtering(byte[] img, int width, int height) {
		Complex[] ft = ft(img, width, height);
		ft[0] = new Complex();
		idft(ft, img, width, height);
	}

}
