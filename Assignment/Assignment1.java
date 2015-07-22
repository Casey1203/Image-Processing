/**
  * Assignment 1
  * You may only edit and submit this file to moodle before the deadline. 
  */
public class Assignment1 {

  /**
    * Task 1
    * Implement the Fast Fourier Transform (FFT) and display the Fourier spectrum.
    * Results should be equivalent to those obtained in the Workshop 3.
    * The implementation details of the FFT can be obtained from section 4.11 of chapter 4. 
    * 6 Marks
    */
    public void fourierTransform(byte[] img, int width, int height) {
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
