import java.text.DecimalFormat;
import java.util.Arrays;

/**
 * A matrix, implemented as a 1D double array. 
 */
public class Matrix {
	/**
	 * number of rows
	 */
	public int rows;
	/**
	 * number of columns
	 */
	public int cols;
	/**
	 * array that stores the matrix content
	 */
	public double[] data; 

	/**
	 * Creates a new matrix with 0 rows and 0 columns.
	 */
	public Matrix() {
		init();
	}
	
	/**
	 * Creates a new matrix with r rows and c columns.
	 * @param r number of rows
	 * @param c number of columns
	 */
	public Matrix(int r, int c) {
		rows = 0;
		cols = 0;
		data = null;

		setSize(r, c);
	}
	
	/**
	 * Creates a new matrix which is a copy of m.
	 * @param m the matrix that will be copied.
	 */
	public Matrix(Matrix m) {
		rows = 0;
		cols = 0;
		data = null;

		setSize(m.rows(), m.cols());
		System.arraycopy(m.data, 0, this.data, 0, m.data.length);
	}

	/**
	 * Creates a new Matrix of size r,c and initializes the values with m.
	 * @param r
	 * @param c
	 * @param m
	 */
	public Matrix(int r, int c, double[] m) {
		rows = 0;
		cols = 0;
		data = null;

		setSize(r, c);
		System.arraycopy(m, 0, this.data, 0, r*c);
	}
	
	/**
	 * Initialize all values to default;
	 * rows =0
	 * cols =0
	 * data = null
	 */
	protected void init() {
		rows = 0;
		cols = 0;
		data = null;
	}
	
	/**
	 * Set the Size of this matrix to r rows, c columns
	 * @param r number of rows
	 * @param c number of columns
	 * @return
	 */
	public void setSize(int r, int c) {
		// check if both Rows and Cols are valid
		if (r < 1 || c < 1) {
			try {
				throw new Exception ("Number of Rows or Cols not valid: "+ r + " - "+c);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		// check if the current data already has correct size
		if (rows * cols == r * c) {
			rows = r;
			cols = c;
			return;
		}

		init();

		data = new double[r * c];
		rows = r;
		cols = c;
	}

	/**
	 * Set the size of the matrix and initialize every element to v
	 * @param r number of rows
	 * @param c number of columns
	 * @param v initial value of all entries
	 * @return
	 */
	public void setSize(int r, int c, double v) {
		setSize(r, c);
		for (int i = 0; i < rows * cols; i++) {
			data[i] = v;
		}
	}

	/**
	 * Set the size of the vector to l (that is a matrix with l rows and 1 column).
	 * @param l
	 * @return
	 */
	public void setSize(int l) {
		setSize(l, 1);
	}
	
	/**
	 * Set every element to zero.
	 */
	public void setZeros() {
		if (isValid()) {
			Arrays.fill(this.data, 0);
		}
	}

	/**
	 * returns the number of rows.
	 * @return number of rows
	 */
	public int rows() {
		return rows;
	}

	/**
	 * return number of columns
	 * @return number of columns
	 */
	public int cols() {
		return cols;
	}

	/** 
	 * return the length of the vector
	 * @return length of the vector
	 */
	public int length() {
		return rows * cols;
	}

	/**
	 * return true if rows !=0;
	 * @return true if rows !=  0
	 */
	public boolean isValid() {
		return rows != 0;
	}
	
	/**
	 * returns a specified row as a double array
	 * @param r row number 
	 * @return the row r as a double array
	 */
	public double[] row(int r) {
		double[] tmp = new double[cols];
		System.arraycopy(data, r * cols, tmp, 0, cols);
		return tmp;
	}

	/**
	 * Retrieves a sub-matrix
	 * @param r
	 * @param rr
	 * @param c
	 * @param cc
	 * @param M
	 * @return
	 */
	public void subMat(int r, int rr, int c, int cc, Matrix M) {
		if (r >= 0 && rr < rows && rr >= r && c >= 0 && cc < cols && cc >= c) {
			M.setSize(rr - r + 1, cc - c + 1);
			for (int i = r; i <= rr; i++) {
				for (int j = c; j <= cc; j++) {
					M.set(i - r, j - c, (this).get(i, j));
				}
			}		
		} else {
			try {
				throw new Exception("Parameter not valid");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	/**
	 * Retrieves a sub-matrix
	 * @param r
	 * @param rr
	 * @param c
	 * @param cc
	 * @return the sub matrix
	 */
	public Matrix subMat(int r, int rr, int c, int cc) {
		Matrix M = new Matrix();
		if (r >= 0 && rr < rows && rr >= r && c >= 0 && cc < cols && cc >= c) {
			M.setSize(rr - r + 1, cc - c + 1);
			for (int i = r; i <= rr; i++) {
				for (int j = c; j <= cc; j++) {
					// M(i-r,j-c) = (*this)(i,j);
					M.set(i - r, j - c, (this).get(i, j));
				}
			}
		}
		return M;
	}

	/**
	 * Retrieves a row as a separate vector
	 * @param r row 
	 * @param M the row vector
	 * @return  
	 */
	public void getRow(int r, Matrix M) {
		subMat(r, r, 0, cols - 1, M);
	}

	/**
	 * Retrieves a row as a separate vector
	 * @param r row 
	 * @return the row vector  
	 */
	public Matrix getRow(int r) {
		return (this).get(r, r, 0, cols - 1);
	}

	/**
	 * Retrieves a column as a separate vector
	 * @param c column
	 * @param M the columns vector
	 * @return
	 */
	public void getCol(int c, Matrix M) {
		subMat(0, rows - 1, c, c, M);
	}

	/**
	 * Retrieves a column as a separate vector
	 * @param c column
	 * @return the column vector.
	 */
	public Matrix getCol(int c) {
		return (this).get(0, rows - 1, c, c);
	}

	/**	
	 * returns the trace (sum of the diagonal elements)
	 * @return the trace (sum of the diagonal elements)
	 */
	public double trace() {
		double sum = (double) 0;

		int l = Math.min(rows, cols);

		for (int i = 0; i < l; i++) {
			sum += (this).get(i, i);
		}

		return sum;
	}

	/**
	 * Performs QR decomposition.
	 * Every MxN matrix A with linearly independent columns can be factored into A=QR; 
	 * The columns of Q are orthonormal, and
	 * The matrix R is upper triangular
	 * @param Q
	 * @param R
	 * @return
	 */
	public void QR(Matrix Q, Matrix R) {
		// check if this matrix has been initialized properly
		if (rows == 0) {
			try {
				throw new Exception("Cannot perfrom QR on an empty matrix");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		// set the sizes of Q and R
		Q.setSize(rows, cols);
		R.setSize(cols, cols);
		
		Q.set(this);
		R.setZeros();

		int i, j, k;

		for (j = 0; j < cols; j++) {
			double norm, val, scale;

			// find the largest element in the current column
			scale = 0.0;
			for (i = 0; i < rows; i++) {
				scale = Math.max(scale, Math.abs((double) (Q.get(i, j))));
			}

			if (scale == 0.0) {
				continue;
			}

			// compute the vector norm for the current column
			norm = 0.0;
			for (i = 0; i < rows; i++) {
				val = Q.get(i, j) / scale;
				norm += val * val;
			}
			norm = Math.sqrt(norm) * scale;

			// assign the norm to the diagonals of R
			// R(j,j) = (Type) norm;
			R.set(j, j, (double) norm);

			// make the current column a unit vector
			for (i = 0; i < rows; i++) {
				// Q(i,j) = (Type) ((Q(i,j))/norm);
				Q.set(i, j, (double) ((Q.get(i, j)) / norm));
			}

			for (k = j + 1; k < cols; k++) {
				// compute the projection onto the current column
				norm = 0.0;
				for (i = 0; i < rows; i++) {
					norm += (double) (Q.get(i, j)) * (double) (Q.get(i, k));
				}

				// update R
				R.set(j, k, (double) norm);

				// remove the component due to the current column
				for (i = 0; i < rows; i++) {
					Q.set(i, k, Q.get(i, k) - (double) (norm * (double) (Q.get(i, j))));
				}
			}
		}
	}

	/**
	 * perform modified QR decomposition such that R is lower-triangular
	 * @param Q
	 * @param R
	 * @return
	 */
	public boolean QR2(Matrix Q, Matrix R) {
		// check if this matrix has been initialized properly
		if (rows == 0) {
			return false;
		}

		// set the sizes of Q and R
		Q.setSize(rows, cols);
		R.setSize(cols, cols);

		// Q = *this;
		Q.set(this);
		R.setZeros();

		int i, j, k;

		for (j = cols - 1; j >= 0; j--) {
			double norm, val, scale;

			// find the largest element in the current column
			scale = 0.0;
			for (i = 0; i < rows; i++) {
				scale = Math.max(scale, Math.abs((double) (Q.get(i, j))));
			}

			if (scale == 0.0) {
				continue;
			}

			// compute the vector norm for the current column
			norm = 0.0;
			for (i = 0; i < rows; i++) {
				val = Q.get(i, j) / scale;
				norm += val * val;
			}
			norm = Math.sqrt(norm) * scale;

			// assign the norm to the diagonals of R
			R.set(j, j, (double) norm);

			// make the current column a unit vector
			for (i = 0; i < rows; i++) {
				Q.set(i, j, (double) ((Q.get(i, j)) / norm));
			}

			for (k = j - 1; k >= 0; k--) {
				// compute the projection onto the current column
				norm = 0.0;
				for (i = 0; i < rows; i++) {
					norm += (double) (Q.get(i, j)) * (double) (Q.get(i, k));
				}

				// update R
				R.set(j, k, (double) norm);

				// remove the component due to the current column
				for (i = 0; i < rows; i++) {
					Q.set(i, k, Q.get(i, k) - (double) (norm * (double) (Q.get(i, j))));
				}
			}
		}

		return true;
	}

	/**
	 * Perform singular value decomposition.
	 * Any MxN real matrix A can be factored by SVD into U=UDV^T.
	 * This is the economy version which will result in
	 * A(MxN)= (MxN) * (Nx1) * (N*N) 
	 * @param U
	 * @param D
	 * @param V
	 * @return
	 */
	public void SVD(Matrix U, Matrix D, Matrix V) {
		int flag, i, its, j, jj, k, l = 0, nm = 0;
		double c, f, h, s, x, y, z;
		double anorm = 0.0, g = 0.0, scale = 0.0;

		int m = Math.max(rows, cols);
		int n = cols;

		// initialize the result matrices
		U.setSize(m, n);
		D.setSize(n);
		V.setSize(n, n);
		U.setZeros();

		// copy the content of this matrix to the upper left submatrix of U
		System.arraycopy((this).data, 0, U.data, 0, rows * cols);

		double[] rv1 = new double[n];

		for (i = 0; i < n; i++) {
			l = i + 1;
			rv1[i] = (double) (scale * g);
			g = s = scale = 0.0;
			if (i < m) {
				for (k = i; k < m; k++)
					scale += Math.abs((double) U.get(k, i));
				if (scale != 0.0) {
					for (k = i; k < m; k++) {
						U.set(k, i, (double) (U.get(k, i) / scale));
						s += (double) (U.get(k, i) * U.get(k, i));
					}
					f = (double) (U.get(i, i));
					g = -Math.copySign(Math.sqrt(s), f);
					h = f * g - s;
					U.set(i, i, (double) (f - g));
					if (i != n - 1) {
						for (j = l; j < n; j++) {
							for (s = 0.0, k = i; k < m; k++)
								s += (double) (U.get(k, i) * U.get(k, j));
							f = s / h;
							for (k = i; k < m; k++)
								U.set(k, j, U.get(k, j) + (double) (f * U.get(k, i)));
						}
					}
					for (k = i; k < m; k++)
						U.set(k, i, (double) (U.get(k, i) * scale));
				}
			}

			D.set(i, (double) (scale * g));
			g = s = scale = 0.0;
			if (i < m && i != n - 1) {
				for (k = l; k < n; k++)
					scale += Math.abs((double) U.get(i, k));
				if (scale != 0.0) {
					for (k = l; k < n; k++) {
						U.set(i, k, (double) (U.get(i, k) / scale));
						s = s + (double) (U.get(i, k) * U.get(i, k));
					}
					f = (double) (U.get(i, l));
					g = -Math.copySign(Math.sqrt(s), f);
					h = f * g - s;
					U.set(i, l, (double) (f - g));
					for (k = l; k < n; k++)
						rv1[k] = (double) (U.get(i, k) / h);
					if (i != m - 1) {
						for (j = l; j < m; j++) {
							for (s = 0.0, k = l; k < n; k++)
								s += (double) (U.get(j, k) * U.get(i, k));
							for (k = l; k < n; k++)
								U.set(j, k, (double) (U.get(j, k) + s * rv1[k]));
						}
					}
					for (k = l; k < n; k++)
						U.set(i, k, (double) (U.get(i, k) * scale));
				}
			}
			anorm = Math.max(anorm, Math.abs((double) D.get(i)) + Math.abs((double) rv1[i]));
		}

		for (i = n - 1; i >= 0; i--) {
			if (i < n - 1) {
				if (g != 0.0) {
					for (j = l; j < n; j++)
						V.set(j, i, (double) ((U.get(i, j) / U.get(i, l)) / g));
					for (j = l; j < n; j++) {
						for (s = 0.0, k = l; k < n; k++)
							s += (double) (U.get(i, k) * V.get(k, j));
						for (k = l; k < n; k++)
							V.set(k, j, (double) (V.get(k, j) + s * V.get(k, i)));
					}
				}
				for (j = l; j < n; j++) {
					// V(i,j) = V(j,i) = (Type) 0.0;
					V.set(i, j, 0.0);
					V.set(j, i, 0.0);
				}

			}
			V.set(i, i, (double) 1.0);
			g = (double) (rv1[i]);
			l = i;
		}
		for (i = n - 1; i >= 0; i--) {
			l = i + 1;
			g = D.get(i);
			if (i < n - 1)
				for (j = l; j < n; j++)
					U.set(i, j, 0.0);
			if (g != 0.0) {
				g = 1.0 / g;
				if (i != n - 1) {
					for (j = l; j < n; j++) {
						for (s = 0.0, k = l; k < m; k++)
							s = s + (double) (U.get(k, i) * U.get(k, j));
						f = (s / U.get(i, i)) * g;
						for (k = i; k < m; k++)
							U.set(k, j, (double) (U.get(k, j) + f * U.get(k, i)));
					}
				}
				for (j = i; j < m; j++)
					U.set(j, i, (double) (U.get(j, i) * g));
			} else {
				for (j = i; j < m; j++)
					U.set(j, i, (double) 0.0);
			}
			// ++U(i,i);
			U.set(i, i, U.get(i, i) + 1);
		}
		for (k = n - 1; k >= 0; k--) {
			for (its = 0; its < 30; its++) {
				flag = 1;
				for (l = k; l >= 0; l--) {
					nm = l - 1;
					if (Math.abs((double) rv1[l]) + anorm == anorm) {
						flag = 0;
						break;
					}
					if (Math.abs((double) D.get(nm)) + anorm == anorm)
						break;
				}
				if (flag != 0) {
					c = 0.0;
					s = 1.0;
					for (i = l; i <= k; i++) {
						f = s * rv1[i];
						if (Math.abs((double) f) + anorm != anorm) {
							g = (double) (D.get(i));
							h = pythag(f, g);
							D.set(i, h);
							h = 1.0 / h;
							c = g * h;
							s = -f * h;
							for (j = 0; j < m; j++) {
								y = (double) (U.get(j, nm));
								z = (double) (U.get(j, i));
								U.set(j, nm, (double) (y * c + z * s));
								U.set(j, i, (double) (z * c - y * s));
							}
						}
					}
				}
				z = D.get(k);
				if (l == k) {
					if (z < 0.0) {
						D.set(k, (double) (-z));
						for (j = 0; j < n; j++)
							V.set(j, k, -V.get(j, k));
					}
					break;
				}

				x = (double) (D.get(l));
				nm = k - 1;
				y = (double) (D.get(nm));
				g = (double) (rv1[nm]);
				h = (double) (rv1[k]);
				f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
				g = pythag(f, 1.0);
				f = ((x - z) * (x + z) + h * ((y / (f + Math.copySign(g, f))) - h)) / x;
				c = s = 1.0;
				for (j = l; j <= nm; j++) {
					i = j + 1;
					g = (double) (rv1[i]);
					y = (double) (D.get(i));
					h = s * g;
					g = c * g;
					z = pythag(f, h);
					rv1[j] = (double) z;
					c = f / z;
					s = h / z;
					f = x * c + g * s;
					g = g * c - x * s;
					h = y * s;
					y = y * c;
					for (jj = 0; jj < n; jj++) {
						x = (double) (V.get(jj, j));
						z = (double) (V.get(jj, i));
						V.set(jj, j, (double) (x * c + z * s));
						V.set(jj, i, (double) (z * c - x * s));
					}
					z = pythag(f, h);
					D.set(j, (double) z);
					if (z != 0.0) {
						z = 1.0 / z;
						c = f * z;
						s = h * z;
					}
					f = (c * g) + (s * y);
					x = (c * y) - (s * g);
					for (jj = 0; jj < m; jj++) {
						y = (double) (U.get(jj, j));
						z = (double) (U.get(jj, i));
						U.set(jj, j, (double) (y * c + z * s));
						U.set(jj, i, (double) (z * c - y * s));
					}
				}
				rv1[l] = (double) 0.0;
				rv1[k] = (double) f;
				D.set(k, (double) x);
			}
		}

		double wmax = 0;
		for (i = 0; i < n; i++) {
			if (D.get(i) > wmax)
				wmax = (double) (D.get(i));
		}

		double singular_value_cutoff = (1e-20);
		double wmin = wmax * singular_value_cutoff;
		for (i = 0; i < n; i++) {
			if (D.get(i) < wmin)
				D.set(i, 0.0);
		}
	}

	/**
	 * Perform singular value decomposition with descending singular values in D
	 * @param U
	 * @param D
	 * @param V
	 * @return
	 */
	public void SVD2(Matrix U, Matrix D, Matrix V) {
		// perform singular value decomposition
		SVD(U, D, V);

		int l = D.length();

		// sorting is not necessary if there are less than 2 values
		if (l < 2) {
			return;
		}
			

		for (int c = 0; c < l - 1; c++) {
			double wmax = D.get(c);
			int cc = 0;

			for (int i = c + 1; i < l; i++) {
				if (D.get(i) > wmax) {
					wmax = D.get(i);
					cc = i;
				}
			}

			if (cc != 0) {
				U.swapCols(c, cc);
				D.swap(c, cc);
				V.swapCols(c, cc);
			}
		}
	}
	
	/**
	 * returns the transpose
	 * @return the transpose
	 */
	public Matrix transpose() {
		Matrix result = new Matrix();
		result.setSize(cols, rows);
			for (int r = 0; r < rows; r++) {
				for (int c = 0; c < cols; c++) {
					result.set(c, r, (this).get(r, c));
				}
			}
		return result;
	}

	/**
	 * return the inverse of the matrix.
	 * @return the inverse of the matrix
	 */
	public Matrix inverse() {
		Matrix U = new Matrix(), D = new Matrix(), V = new Matrix(), W = new Matrix();

		SVD(U, D, V);

		int l = D.length();

		W.setSize(l, l);

		W.setZeros();

		for (int i = 0; i < l; i++) {
			double d = D.get(i);
			double singular_value_cutoff = (1e-20);
			if (d > singular_value_cutoff) {
				W.set(i, i, (double) (1 / d));
			}
		}

		return V.mul(W.mul(U.transpose()));
	}

	/**
	 * set values of the submatrix to M,
	 * 
	 * @param r
	 * @param c
	 * @param M
	 * @return
	 */
	public boolean set(int r, int c, Matrix M) {
		int rr = r + M.rows() - 1;
		int cc = c + M.cols() - 1;

		if (r >= 0 && rr >= r && rr < rows && c >= 0 && cc >= c && cc < cols) {
			for (int i = r; i <= rr; i++) {
				for (int j = c; j <= cc; j++) {
					(this).set(i, j, M.get(i - r, j - c));
				}
			}

			return true;
		} else {
			return false;
		}
	}

	/**
	 * set values of the submatrix to M(r1:r2,c1:c2)
	 * @param r
	 * @param c
	 * @param M
	 * @param r1
	 * @param r2
	 * @param c1
	 * @param c2
	 * @return
	 */
	public boolean set(int r, int c, Matrix M, int r1, int r2, int c1, int c2) {
		int rr = r + r2 - r1;
		int cc = c + c2 - c1;

		if (r1 >= 0 && r2 >= r1 && r2 < M.rows() && c1 >= 0 && c2 >= c1 && c2 < M.cols() && r >= 0 && rr >= r && rr < rows
				&& c >= 0 && cc >= c && cc < cols) {
			for (int i = r1; i <= r2; i++) {
				for (int j = c1; j <= c2; j++) {
					(this).set(r + i - r1, c + j - c1, M.get(i, j));
				}
			}

			return true;
		} else {
			return false;
		}
	}

	/**
	 * 	// set values of the submatrix to T,
	 * @param r
	 * @param rr
	 * @param c
	 * @param cc
	 * @param T
	 * @return
	 */
	public boolean set(int r, int rr, int c, int cc, double T) {
		if (r >= 0 && rr >= r && rr < rows && c >= 0 && cc >= c && cc < cols) {
			for (int i = r; i <= rr; i++) {
				for (int j = c; j <= cc; j++) {
					(this).set(i, j, T);
				}
			}

			return true;
		} else {
			return false;
		}
	}

	/**
	 * 	// set the row of M,
	 * @param r
	 * @param Row
	 * @return
	 */
	public boolean setRow(int r, Matrix Row) {
		if (r >= 0 && r < rows && Row.rows() == 1 && Row.cols() == cols) {
			for (int c = 0; c < cols; c++) {
				(this).set(r, c, Row.get(c));
			}

			return true;
		} else {
			return false;
		}
	}

	/**
	 * 	set the column of M
	 * @param c
	 * @param Col
	 * @return
	 */
	public boolean setCol(int c, Matrix Col) {
		if (c >= 0 && c < cols && Col.cols() == 1 && Col.rows() == rows) {
			for (int r = 0; r < rows; r++) {
				(this).set(r, c, Col.get(r));
			}

			return true;
		} else {
			return false;
		}
	}

	/**
	 * swap rows r and rr
	 * @param r
	 * @param rr
	 * @return
	 */
	public boolean swapRows(int r, int rr) {
		if (r == rr) {
			return true;
		} else if (r >= 0 && r < rows && rr >= 0 && rr < rows) {
			double temp;
			for (int c = 0; c < cols; c++) {
				temp = (this).get(r, c);
				(this).set(r, c, (this).get(rr, c));
				(this).set(rr, c, temp);
			}

			return true;
		} else {
			return false;
		}
	}

	/**
	 * 	swap columns c and cc
	 * @param c
	 * @param cc
	 * @return
	 */
	public boolean swapCols(int c, int cc) {
		if (c == cc) {
			return true;
		} else if (c >= 0 && c < cols && cc >= 0 && cc < cols) {
			double temp;
			for (int r = 0; r < rows; r++) {
				temp = (this).get(r, c);
				(this).set(r, c, (this).get(r, cc));
				(this).set(r, cc, temp);
			}

			return true;
		} else {
			return false;
		}
	}

	/**
	 * swap values at location i and j
	 * @param i
	 * @param j
	 * @return
	 */
	public boolean swap(int i, int j) {
		if (i == j) {
			return true;
		} else if (i >= 0 && i < length() && j >= 0 && j < length()) {
			double temp;

			temp = data[i];
			data[i] = data[j];
			data[j] = temp;

			return true;
		} else {
			return false;
		}
	}

		
	/**
	 * Gets element at position (r,c)
	 * @param r row
	 * @param c column
	 * @return
	 */
	public double get(int r, int c) {
		return data[r * cols + c];
	}

	/**
	 * Sets the element at position (r,c) to v
	 * @param r
	 * @param c
	 * @param v
	 */
	public void set(int r, int c, double v) {
		data[r * cols + c] = v;
	}

	/**
	 * 	access M(i) as a vector without range checking
	 * @param i
	 * @return
	 */
	public double get(int i) {
		return data[i];
	}
	
	/**
	 * set M(i) = v;
	 * @param i
	 * @param v
	 */
	public void set(int i, double v) {
		data[i] = v;
	}

	/**
	 * 	return the submatrix (r:rr,c:cc)
	 * @param r
	 * @param rr
	 * @param c
	 * @param cc
	 * @return
	 */
	public Matrix get(int r, int rr, int c, int cc) {
		Matrix result = new Matrix();
		if (r >= 0 && rr >= r && rr < rows && c >= 0 && cc >= c && cc < cols) {
			result.setSize(rr - r + 1, cc - c + 1);
			for (int i = r; i <= rr; i++) {
				for (int j = c; j <= cc; j++) {
					result.set(i - r, j - c, (this).get(i, j));
				}
			}
		} else {
			try {
				throw new Exception("not valid parameter");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return result;
	}

	/**
	 * +
	 * 
	 * @param M
	 * @return
	 */
	public Matrix plus(Matrix M) {
		Matrix result = new Matrix();
		if (rows == M.rows() && cols == M.cols()) {
			result.setSize(rows, cols);
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					result.set(i, j, (this).get(i, j) + M.get(i, j));
				}
			}
		} else {
			try {
				throw new Exception("not valid parameter");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return result;
	}
	/**
	 * -
	 * @param M
	 * @return
	 */
	public Matrix minus(Matrix M) {
		Matrix result = new Matrix();
		if (rows == M.rows() && cols == M.cols()) {
			result.setSize(rows, cols);
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					result.set(i, j, get(i, j) - M.get(i, j));
				}
			}
		} else {
			try {
				throw new Exception("not valid parameter");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return result;
	}

	/**
	 * *
	 * @param M
	 * @return
	 */
	public Matrix mul(Matrix M) {
		Matrix result = null;
		if (cols == M.rows()) {
			result = new Matrix(rows, M.cols());
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < M.cols(); j++) {
					for (int k = 0; k < cols; k++) {
						result.set(i, j, (result.get(i, j) + (this).get(i, k) * M.get(k, j)));
					}
				}
			}
		} else {
			try {
				throw new Exception("not valid parameter: "+this.cols+" and "+M.rows());
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return result;
	}

	/**
	 * +
	 * 
	 * @param T
	 * @return
	 */
	public Matrix plus(double T) {
		Matrix result = new Matrix();
		result.setSize(rows, cols);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				result.set(i, j, get(i, j) + T);
			}
		}
		return result;
	}

	/**
	 * -
	 * 
	 * @param T
	 * @return
	 */
	public Matrix minus(double T) {
		Matrix result = new Matrix();
		result.setSize(rows, cols);
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					result.set(i, j, get(i, j) - T);
				}
			}
		return result;
	}

	/**
	 * *
	 * @param T
	 * @return
	 */
	public Matrix mul(double T) {
		Matrix result = new Matrix();
		result.setSize(rows, cols);
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					result.set(i, j, get(i, j) * T);
				}
			}
		return result;
	}

	/**
	 * \
	 * @param T
	 * @return
	 */
	public Matrix div(double T) {
		Matrix result = new Matrix();
		result.setSize(rows, cols);
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					result.set(i, j, get(i, j) / T);
				}
			}
		return result;
	}


	/**
	 * Assign matrix M to this matrix
	 * @param M
	 * @return
	 */
	public Matrix set(Matrix M) {
		if (this != M) {
			setSize(M.rows(), M.cols());
			System.arraycopy(M.data, 0, this.data, 0, M.data.length);
		}
		return this;
	}


	private double pythag(double a, double b) {
		double at, bt;

		at = Math.abs(a);
		bt = Math.abs(b);

		if (at > bt)
			return (at * Math.sqrt(1.0 + (bt / at) * (bt / at)));
		else if (bt == 0.0)
			return 0.0;
		else
			return (bt * Math.sqrt(1.0 + (at / bt) * (at / bt)));
	}


	/**
	 * Returns a string representation of the matrix.
	 */
	public String toString() {
		DecimalFormat df2 = new DecimalFormat("+0000.00;-0000.00");
		String result = new String("");
		for (int r = 0; r < rows; r++) {
			result = result.concat("|");
			for (int c = 0; c < cols; c++) {
				result = result.concat(df2.format(data[r * cols + c]) + " ");
			}
			if (r==rows-1) {
				result = result.concat("|");	
			} else {
				result = result.concat("|" + "\n");	
			}
			
		}
		return result;
	}
}
