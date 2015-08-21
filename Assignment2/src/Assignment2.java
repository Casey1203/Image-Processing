import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;

/**
  * Assignment 2
  * You may only edit and submit this file to moodle before the deadline. 
  */
public class Assignment2 {
   /**
    * Task 1
    * Implement 2D Gaussian smoothing. Hints:
    * - Perform 2D filtering by applying 1D filters. 
    * - Compute and use a proper filter size for a 1D Gaussian mask based on the sigma parameter. 
    * - Employ partial filters to handle the boundary of the image 
    * 5 Marks
    */
	public void gaussianSmooth(final byte[] f, int width, int height, double sigma) {
		double[] new_f = new double[f.length];
		for(int i = 0; i < f.length; i++){
			new_f[i] = f[i] & 0xff;
		}
		int filterSize = 6 * (int)sigma + 1;//has relation with sigma
		int halfFilterSize = filterSize / 2;
		double[] filter = new double[filterSize];
		double sumOfFilter = 0.0;
		//1D gaussian filter
		for(int x = 0; x < filterSize; x++){
			filter[x] = 1/(Math.sqrt(2.0 * Math.PI) * sigma) * Math.exp(-Math.pow((x - halfFilterSize)*1.0/sigma, 2)/2.0);
			sumOfFilter += filter[x];
		}
		//scaling
		for(int i = 0; i < filter.length; i++){
			filter[i] = filter[i] / sumOfFilter;
		}
		convolve(new_f, filter, filterSize, width, height);
		for(int i = 0; i < new_f.length; i++){
			f[i] = (byte)new_f[i];
		}
	}

	/**
	 *
	 * @param f double type
	 * @param width
	 * @param height
	 * @param sigma
	 */
	private void gaussianSmoothdoubleInput(double[] f, int width, int height, double sigma){
		int filterSize = (int)sigma * 6 + 1;
		int halfFilterSize = filterSize / 2;
		double[] filter = new double[filterSize];
		double sumOfFilter = 0.0;
		//1D gaussian filter
		for(int x = 0; x < filterSize; x++){
			filter[x] = 1/(Math.sqrt(2.0 * Math.PI) * sigma) * Math.exp(-Math.pow((x - halfFilterSize)*1.0/sigma, 2)/2.0);
			sumOfFilter += filter[x];
		}
		//scaling
		for(int i = 0; i < filter.length; i++){
			filter[i] = filter[i] / sumOfFilter;
		}
		convolve(f, filter, filterSize, width, height);
	}
	private void convolve(double[] img, double[] filter, int filterSize, int width, int height){
		double[] new_img = img.clone();
		int halfFilterSize = filterSize / 2;
		for(int x = halfFilterSize; x < height-halfFilterSize; x++){
			for(int y = halfFilterSize; y < width-halfFilterSize; y++){
				img[x * width + y] = 0;
				for(int a = -halfFilterSize; a <= halfFilterSize; a++){
					double innerSum = 0;
					for(int b = -halfFilterSize; b <= halfFilterSize; b++){
						innerSum += new_img[(x - a) * width + y - b] * filter[b + halfFilterSize];
					}
					img[x * width + y] += filter[a + halfFilterSize] * innerSum;
				}
			}
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
	 * Implement the Harris Corner Detection Algorithm. Hints:
	 * - Make proper use of the sigma and threshold values.
	 * - Suppress the non-maximal corner candidates
	 * - Use Gaussian smoothing on the squared derivative images
	 * - Computer corners up to sub-pixel accuracy
	 * 5 Marks
	 */
	public ArrayList<double[]> detectCorners(final byte[] f, int width, int height, double sigma, double threshold) {
		double[] fx2 = new double[f.length];
		double[] fy2 = new double[f.length];
		double[] fxy = new double[f.length];
		double[] R = new double[f.length];

		for(int x = 1; x < height-1; x++){
			for(int y = 1; y < width-1; y++){
				int fx = ((int)f[(x+1) * width + y]&0xff) - ((int)f[(x-1) * width + y] & 0xff);
				int fy = ((int)f[x * width + y + 1]&0xff) - ((int)f[x * width + y - 1] & 0xff);
				fx2[x * width + y] = fx * fx;
				fy2[x * width + y] = fy * fy;
				fxy[x * width + y] = fx * fy;
			}
		}
		gaussianSmoothdoubleInput(fx2, width, height, sigma);
		gaussianSmoothdoubleInput(fy2, width, height, sigma);
		gaussianSmoothdoubleInput(fxy, width, height, sigma);

		// calculate R
		for(int i = 0; i < R.length; i ++){
			R[i] = (fx2[i] * fy2[i] - fxy[i] * fxy[i]) - 0.04 * ((fx2[i] + fy2[i]) * (fx2[i] + fy2[i]));
		}

		ArrayList<double[]> cornersOut = new ArrayList<double[]>();

		int maskSize = 3;
		int halfMaskSize = maskSize/2;
		for(int x = halfMaskSize; x < height-halfMaskSize; x++){
			for(int y = halfMaskSize; y < width-halfMaskSize; y++){
				if(isLocalMaxima(R, x, y, width, height, maskSize) && R[x * width + y] > threshold){
					//parameter defined by myself
					double a = (R[(x+1) * width+y] + R[(x-1) * width+y] - 2 * R[x * width + y])/2.0;
					double b = (R[x * width + y + 1] + R[x * width + y - 1] - 2 * R[x * width + y])/2.0;
					double c = (R[(x+1) * width + y + 1] + R[(x-1) * width + y - 1])/2.0 - R[x * width + y] - a - b;
					double d = R[(x+1) * width + y] - R[x * width + y] - a * (2 * x + 1) - c * y;
					double e = R[x * width + y + 1] - R[x * width + y - 1] - b * (2 * y + 1) - c * x;

					double x_appr = (c * e - 2 * b * d) / (4 * a * b - c * c);
					double y_appr = (c * d - 2 * a * e) / (4 * a * b - c * c);
					double[] sub_acc = new double[2];
					sub_acc[0] = y_appr;
					sub_acc[1] = x_appr;
					cornersOut.add(sub_acc);
				}
			}
		}

		return cornersOut;
	}
	private ArrayList<Double> getNeighborNbyN(double[] img, int x, int y, int width, int height, int filterSize){
		int halfFilterSize = (filterSize-1)/2;
		ArrayList<Double> effAdjIntensity = new ArrayList<Double>();
		for(int i = -halfFilterSize; i <= halfFilterSize; i++){
			for(int j = -halfFilterSize; j <= halfFilterSize; j++){
				int x_neighbor = x+i;
				int y_neighbor = y+j;
				if(i != 0 && j != 0){//remove itself
					effAdjIntensity.add(img[x_neighbor * width + y_neighbor]);
				}
			}
		}
		return effAdjIntensity;
	}
	private boolean isLocalMaxima(double[] img, int x, int y, int width, int height, int maskSize){
		ArrayList<Double> effAdjIntensity = getNeighborNbyN(img, x, y, width, height, maskSize);
		double max = Collections.max(effAdjIntensity);
		if(img[x * width + y] > max){
			return true;
		}
		else{
			return false;
		}
	}
	/**
	  * This method is called when the Task 3 is clicked in the UI.
	  * This method together with various helper methods will generate many 2D to 3D correspondences. 
	  * No need to edit this function. 
	  */
	public Matrix determineProjectionMatrix(ArrayList<double[]> corners, ArrayList<double[]> pointList2D, ArrayList<double[]> pointList3D) {
		Matrix projectionMatrix = new Matrix();
		ArrayList<double[]> yzPlane2DPoints = new ArrayList<double[]>();
		ArrayList<double[]> yzPlane3DPoints = new ArrayList<double[]>();
		ArrayList<double[]> xzPlane2DPoints = new ArrayList<double[]>();
		ArrayList<double[]> xzPlane3DPoints = new ArrayList<double[]>();
		for (int i=0;i<pointList3D.size();i++) {
			if (pointList3D.get(i)[1] == 0) {
				xzPlane2DPoints.add(pointList2D.get(i));
				xzPlane3DPoints.add(pointList3D.get(i));
			} else if (pointList3D.get(i)[0] == 0) {
				yzPlane2DPoints.add(pointList2D.get(i));
				yzPlane3DPoints.add(pointList3D.get(i));
			}
		}
		ArrayList<double[]> final2D = new ArrayList<double[]>();
		ArrayList<double[]> final3D = new ArrayList<double[]>();
		if ((xzPlane2DPoints.size() >= 4 && yzPlane2DPoints.size() >= 4)) {
			double gridSize = 0.5;
			double max_u, max_v;
			int[] gridCount = new int[3];
			if (yzPlane2DPoints.size() >= 4) {
				Matrix P_YZ = performPlanarCalibration(yzPlane2DPoints, yzPlane3DPoints, 0);
				max_u = 0;
				max_v = 0;
				for (int i = 0; i < yzPlane3DPoints.size(); i++) {
					max_u = Math.max(yzPlane3DPoints.get(i)[1], max_u);
					max_v = Math.max(yzPlane3DPoints.get(i)[2], max_v);
				}
				gridCount[0] = 1;
				gridCount[1] = (int) Math.floor(max_u / gridSize + 1.5);
				gridCount[2] = (int) Math.floor(max_v / gridSize + 1.5);
				obtainRefiningPoints(corners, P_YZ, gridSize, gridCount, final2D, final3D);
			}
			if (xzPlane2DPoints.size() >= 4) {
				Matrix P_XZ = performPlanarCalibration(xzPlane2DPoints, xzPlane3DPoints, 1);
				max_u = 0;
				max_v = 0;
				for (int i = 0; i < (int) xzPlane3DPoints.size(); i++) {
					max_u = Math.max(xzPlane3DPoints.get(i)[0], max_u);
					max_v = Math.max(xzPlane3DPoints.get(i)[2], max_v);
				}
				gridCount[0] = (int) Math.floor(max_u / gridSize + 1.5);
				gridCount[1] = 1;
				gridCount[2] = (int) Math.floor(max_v / gridSize + 1.5);
				obtainRefiningPoints(corners, P_XZ, gridSize, gridCount, final2D, final3D);
			}
		}
		corners.clear();
		corners.addAll(final2D);
		pointList2D.clear();
		pointList2D.addAll(final2D);
		pointList3D.clear();
		pointList3D.addAll(final3D);
		return performCalibration(final2D, final3D);
	}

	/**
	  * No need to edit this function. 
	  */
	private Matrix performPlanarCalibration(ArrayList<double[]> points2d_in, ArrayList<double[]> points3d_in, int planeID) {
		int numberOfPoints = points2d_in.size();
		Matrix A = new Matrix(numberOfPoints * 2, 8);
		Matrix B = new Matrix(numberOfPoints * 2, 1);
		int coordinateIndex = 0;
		if (planeID == 0) coordinateIndex = 1;
		if (planeID == 1) coordinateIndex = 0;
		for (int i = 0; i < numberOfPoints; i++) {
			int c = 0;
			A.set(2*i, c++, points3d_in.get(i)[coordinateIndex]);
			A.set(2*i, c++, points3d_in.get(i)[2]);
			A.set(2*i, c++, 1.0);
			A.set(2*i, c++, 0.0);
			A.set(2*i, c++, 0.0);
			A.set(2*i, c++, 0.0);
			A.set(i << 1, c++, -points2d_in.get(i)[0] * points3d_in.get(i)[coordinateIndex]);
			A.set(i << 1, c++, -points2d_in.get(i)[0] * points3d_in.get(i)[2]);
			c = 0;
			A.set(i * 2 + 1, c++, 0.0);
			A.set(i * 2 + 1, c++, 0.0);
			A.set(i * 2 + 1, c++, 0.0);
			A.set(i * 2 + 1, c++, points3d_in.get(i)[coordinateIndex]);
			A.set(i * 2 + 1, c++, points3d_in.get(i)[2]);
			A.set(i * 2 + 1, c++, 1.0);
			A.set(i * 2 + 1, c++, -(points2d_in).get(i)[1] * (points3d_in).get(i)[coordinateIndex]);
			A.set(i * 2 + 1, c++, -(points2d_in).get(i)[1] * (points3d_in).get(i)[2]);
			B.set(i * 2, points2d_in.get(i)[0]);
			B.set(i * 2 + 1, points2d_in.get(i)[1]);
		}
		Matrix x = A.inverse().mul(B);
		Matrix prj = new Matrix(3, 3);
		for (int i = 0; i < 8; i++) {
			prj.set(i / 3, i % 3, x.get(i, 0));
		}
		prj.set(2, 2, 1.0);
		Matrix P = new Matrix(3, 4);
		for (int i = 0; i < 3; i++) {
			int c = 0;
			for (int j = 0; j < 4; j++) {
				if (j == planeID)
					P.set(i, j, 0.0);
				else
					P.set(i, j, prj.get(i, c++));
			}
		}
		return P;
	}
	/**
	  * No need to edit this function.
	  */
	private void obtainRefiningPoints(ArrayList<double[]> points2d_in, Matrix P, double gridSize_in, int gridCount_in[], ArrayList<double[]> points2d_out, ArrayList<double[]> points3d_out) {
		double fuzziness = 25;
		for (int x = (gridCount_in[0] > 1 ? 1 : 0); x < gridCount_in[0]; x += 2) {
			for (int y = (gridCount_in[1] > 1 ? 1 : 0); y < gridCount_in[1]; y += 2) {
				for (int z = 1; z < gridCount_in[2]; z += 2) {
					Matrix x3D = new Matrix(4, 1);
					x3D.set(0, 0, x * gridSize_in);
					x3D.set(1, 0, y * gridSize_in);
					x3D.set(2, 0, z * gridSize_in);
					x3D.set(3, 0, 1.0);
					Matrix x2D = P.mul(x3D);
					x2D = x2D.div(x2D.get(2));
					double minDistance = 100000;
					int minDistanceIndex = 0;
					for (int j=0;j<points2d_in.size();j++) {
						double distance = Math.pow(x2D.get(0) - points2d_in.get(j)[0], 2) + Math.pow(x2D.get(1) - points2d_in.get(j)[1], 2);
						if (distance < minDistance) {
							minDistance = distance;
							minDistanceIndex = j;
						}
					}
					if (minDistance < fuzziness) {
						points2d_out.add(points2d_in.get(minDistanceIndex));
						points3d_out.add(new Matrix(x3D).data);
					}
				}
			}
		}
	}

    /**
      * Task 3:
      * Perform camera calibration and return the 3x4 camera projection matrix from the provided
      * correspondences. Hints:
      * - Solve the equation Ap = 0 on slide 5 of chapter 7. 
      * - The matrix A is of size (points2d.size() * 2, 12)
      * - Use the provided function SVD2 in the class Matrix.
      * - The solution to Ap = 0 will be the last column in V return by SVD2.
      * 5 Marks
      */
	private Matrix performCalibration(ArrayList<double[]> points2d, ArrayList<double[]> points3d) {
		double[] matrixElement = new double[points2d.size() * 2 * 12];
		int matrixwidth = 12;
		for(int i = 0; i < points2d.size(); i++){
			double[] cur3DPoint = points3d.get(i);
			double[] cur2DPoint = points2d.get(i);
			double[] row1 = {cur3DPoint[0], cur3DPoint[1], cur3DPoint[2], 1, 0, 0, 0, 0,
					-cur2DPoint[0] * cur3DPoint[0], -cur2DPoint[0] * cur3DPoint[1], -cur2DPoint[0] * cur3DPoint[2], -cur2DPoint[0]};
			double[] row2 = {0, 0, 0, 0, cur3DPoint[0], cur3DPoint[1], cur3DPoint[2], 1,
					-cur2DPoint[1] * cur3DPoint[0], -cur2DPoint[1] * cur3DPoint[1], -cur2DPoint[1] * cur3DPoint[2], -cur2DPoint[1]};
			for(int j = 0; j < row1.length; j++){
				matrixElement[2 * i * matrixwidth + j] = row1[j];
				matrixElement[(2 * i + 1) * matrixwidth + j] = row2[j];
			}
		}
		Matrix A = new Matrix(points2d.size() * 2, 12, matrixElement);
//		System.out.println(A.toString());
		int rows = A.rows;
		int cols = A.cols;

		Matrix U = new Matrix(rows, rows);
		Matrix V = new Matrix(cols, cols);
		Matrix D = new Matrix(rows, cols);
		A.SVD2(U, D, V);
		Matrix p_col = V.getCol(cols-1);
		double [] p_ele = new double[cols];
		for(int i = 0; i < cols; i++){
			p_ele[i] = p_col.get(i);
			System.out.println(p_ele[i]);
		}
		//scale p
		for(int i = 0; i < p_ele.length; i++){
			p_ele[i] = p_ele[i] / p_ele[p_ele.length-1];
		}
		Matrix P = new Matrix(3, 4, p_ele);
		return P;
	}
}
