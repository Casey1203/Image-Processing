import javax.swing.JOptionPane;


public class Workshop5 {
	/**
	 * Decomposes the projection matrix into the intrinsic parameters (K) and the extrinsic parameters (RT).
	 * @param projectionMatrix_in the projection matrix to be decomposed
	 * @param K_out the camera intrinsic parameters
	 * @param RT_out the camera extrinsic parameters 
	 */
	public static void decomposeProjectionMatrix(Matrix P, final Matrix K, final Matrix R, final Matrix T) {
		Matrix Q = new Matrix(3,3);
		Matrix RPrime = new Matrix(3,3);
		P.subMat(0, 2, 0, 2).inverse().QR(Q, RPrime);

	}


	public static void main(String args[]) {
		Matrix groundTruthK = new Matrix(3, 3, new double[]{
			100, 0, 50,
			0, 200, 60,
			0, 0, 1
		});
		Matrix groundTruthRT = new Matrix (3, 4, new double[]{
			1, 0, 0, 10, 
			0, 1, 0, 20,
			0, 0, 1, 30
		});
		Matrix P = groundTruthK.mul(groundTruthRT);
		Matrix K = new Matrix(3, 3);
		Matrix R = new Matrix(3, 3);
		Matrix T = new Matrix(3, 1);
		Workshop5.decomposeProjectionMatrix(P, K, R, T);
		JOptionPane.showMessageDialog(null, "P = \n" + P + "\n\nK = \n" + K + "\n\nR = \n" + R+ "\n\nT = \n" + T, "Calibration Result", JOptionPane.INFORMATION_MESSAGE);
	}


}
