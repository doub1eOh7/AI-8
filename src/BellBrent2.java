// The contents of this file are dedicated to the public domain.
// (See http://creativecommons.org/publicdomain/zero/1.0/)
	import java.util.Random;
	import java.util.ArrayList;
	import java.util.Iterator;
	import java.awt.image.BufferedImage;
	import java.awt.Color;
	import java.io.File;
	import javax.imageio.ImageIO;
	import java.util.HashMap;
	import java.util.Scanner;
	import java.io.FileNotFoundException;
	import java.io.PrintWriter;
	import java.util.Map;
	import java.lang.StringBuilder;
	
class BellBrent2 implements IAgent
{
	
	// ----------------------------------------------------------------
	// The contents of this file are distributed under the CC0 license.
	// See http://creativecommons.org/publicdomain/zero/1.0/
	// ----------------------------------------------------------------


	/// Provides several useful static methods for operating on arrays of doubles
	public static class Vec
	{
		public static String toString(double[] vec) {
			StringBuilder sb = new StringBuilder();
			sb.append("[");
			if(vec.length > 0) {
				sb.append(Double.toString(vec[0]));
				for(int i = 1; i < vec.length; i++) {
					sb.append(",");
					sb.append(Double.toString(vec[i]));
				}
			}
			sb.append("]");
			return sb.toString();
		}

		public static void setAll(double[] vec, double val) {
			for(int i = 0; i < vec.length; i++)
				vec[i] = val;
		}

		public static double squaredMagnitude(double[] vec) {
			double d = 0.0;
			for(int i = 0; i < vec.length; i++)
				d += vec[i] * vec[i];
			return d;
		}

		public static void normalize(double[] vec) {
			double mag = squaredMagnitude(vec);
			if(mag <= 0.0) {
				setAll(vec, 0.0);
				vec[0] = 1.0;
			} else {
				double s = 1.0 / Math.sqrt(mag);
				for(int i = 0; i < vec.length; i++)
					vec[i] *= s;
			}
		}

		public static void copy(double[] dest, double[] src) {
			if(dest.length != src.length)
				throw new IllegalArgumentException("mismatching sizes");
			for(int i = 0; i < src.length; i++) {
				dest[i] = src[i];
			}
		}

		public static double[] copy(double[] src) {
			double[] dest = new double[src.length];
			for(int i = 0; i < src.length; i++) {
				dest[i] = src[i];
			}
			return dest;
		}

		public static void add(double[] dest, double[] src) {
			if(dest.length != src.length)
				throw new IllegalArgumentException("mismatching sizes");
			for(int i = 0; i < dest.length; i++) {
				dest[i] += src[i];
			}
		}

		public static void scale(double[] dest, double scalar) {
			for(int i = 0; i < dest.length; i++) {
				dest[i] *= scalar;
			}
		}

		public static double dotProduct(double[] a, double[] b) {
			if(a.length != b.length)
				throw new IllegalArgumentException("mismatching sizes");
			double d = 0.0;
			for(int i = 0; i < a.length; i++)
				d += a[i] * b[i];
			return d;
		}

		public static double squaredDistance(double[] a, double[] b) {
			if(a.length != b.length)
				throw new IllegalArgumentException("mismatching sizes");
			double d = 0.0;
			for(int i = 0; i < a.length; i++) {
				double t = a[i] - b[i];
				d += t * t;
			}
			return d;
		}

		public static void clip(double[] vec, double min, double max) {
			if(max < min)
				throw new IllegalArgumentException("max must be >= min");
			for(int i = 0; i < vec.length; i++) {
				vec[i] = Math.max(min, Math.min(max, vec[i]));
			}
		}

		public static double[] concatenate(double[] a, double[] b) {
			double[] c = new double[a.length + b.length];
			for(int i = 0; i < a.length; i++)
				c[i] = a[i];
			for(int i = 0; i < b.length; i++)
				c[a.length + i] = b[i];
			return c;
		}

	}


	/// This stores a matrix, A.K.A. data set, A.K.A. table. Each element is
	/// represented as a double value. Nominal values are represented using their
	/// corresponding zero-indexed enumeration value. For convenience,
	/// the matrix also stores some meta-data which describes the columns (or attributes)
	/// in the matrix.
	public static class Matrix
	{
		/// Used to represent elements in the matrix for which the value is not known.
		public static final double UNKNOWN_VALUE = -1e308; 

		// Data
		private ArrayList<double[]> m_data = new ArrayList<double[]>(); //matrix elements

		// Meta-data
		private String m_filename;                          // the name of the file
		private ArrayList<String> m_attr_name;                 // the name of each attribute (or column)
		private ArrayList<HashMap<String, Integer>> m_str_to_enum; // value to enumeration
		private ArrayList<HashMap<Integer, String>> m_enum_to_str; // enumeration to value

		/// Creates a 0x0 matrix. (Next, to give this matrix some dimensions, you should call:
		///    loadARFF
		///    setSize
		///    addColumn, or
		///    copyMetaData
		@SuppressWarnings("unchecked")
		public Matrix() 
		{
			this.m_filename    = "";
			this.m_attr_name   = new ArrayList<String>();
			this.m_str_to_enum = new ArrayList<HashMap<String, Integer>>();
			this.m_enum_to_str = new ArrayList<HashMap<Integer, String>>();
		}

		public Matrix(int rows, int cols)
		{
			this.m_filename    = "";
			this.m_attr_name   = new ArrayList<String>();
			this.m_str_to_enum = new ArrayList<HashMap<String, Integer>>();
			this.m_enum_to_str = new ArrayList<HashMap<Integer, String>>();
			setSize(rows, cols);
		}

		public Matrix(Matrix that)
		{
			setSize(that.rows(), that.cols());
			m_filename = that.m_filename;
			m_attr_name = new ArrayList<String>();
			m_str_to_enum = new ArrayList<HashMap<String, Integer>>();
			m_enum_to_str = new ArrayList<HashMap<Integer, String>>();
			copyBlock(0, 0, that, 0, 0, that.rows(), that.cols());
		}

		/// Loads the matrix from an ARFF file
		public void loadARFF(String filename)
		{
			HashMap<String, Integer> tempMap  = new HashMap<String, Integer>(); //temp map for int->string map (attrInts)
			HashMap<Integer, String> tempMapS = new HashMap<Integer, String>(); //temp map for string->int map (attrString)
			int attrCount = 0; // Count number of attributes
			int lineNum = 0; // Used for exception messages
			Scanner s = null;
			m_str_to_enum.clear();
			m_enum_to_str.clear();
			m_attr_name.clear();

			try
			{
				s = new Scanner(new File(filename));
				while (s.hasNextLine())
				{
					lineNum++;
					String line  = s.nextLine().trim();
					String upper = line.toUpperCase();

					if (upper.startsWith("@RELATION"))
						m_filename = line.split(" ")[1];
					else if (upper.startsWith("@ATTRIBUTE"))
					{
						String[] pieces = line.split("\\s+");
						m_attr_name.add(pieces[1]);
						
						tempMap.clear();
						tempMapS.clear();
						
						// If the attribute is nominal
						if (pieces[2].startsWith("{"))
						{
							// Splits this string based on curly brackets or commas
							String[] attributeNames = pieces[2].split("[{},]");
							int valCount = 0;
							
							for (String attribute : attributeNames)
							{
								if (!attribute.equals("")) // Ignore empty strings
								{
									tempMapS.put(valCount, attribute);
									tempMap.put(attribute, valCount++);
								}
							}
						}
						
						// The attribute is continuous if it wasn't picked up in the previous "if" statement
						
						m_str_to_enum.add(new HashMap<String, Integer>(tempMap));
						m_enum_to_str.add(new HashMap<Integer, String>(tempMapS));
						
						attrCount++;
					}
					else if (upper.startsWith("@DATA"))
					{
						m_data.clear();
						
						while (s.hasNextLine())
						{
							double[] temp = new double[attrCount];

							lineNum++;
							line  = s.nextLine().trim();
							
							if (line.startsWith("%") || line.isEmpty()) continue;
							String[] pieces = line.split(",");
							
							if (pieces.length < attrCount) throw new IllegalArgumentException("Expected more elements on line: " + lineNum + ".");
							
							for (int i = 0; i < attrCount; i++)
							{
								int vals   = valueCount(i);
								String val = pieces[i];
								
								// Unknown values are always set to UNKNOWN_VALUE
								if (val.equals("?"))
								{
									temp[i] = UNKNOWN_VALUE;
									continue;
								}
			
								// If the attribute is nominal
								if (vals > 0)
								{
									HashMap<String, Integer> enumMap = m_str_to_enum.get(i);
									if (!enumMap.containsKey(val))
										throw new IllegalArgumentException("Unrecognized enumeration value " + val + " on line: " + lineNum + ".");
										
									temp[i] = (double)enumMap.get(val);
								}
								else
									temp[i] = Double.parseDouble(val); // The attribute is continuous
							}
							
							m_data.add(temp);
						}
					}
				}
			}
			catch (FileNotFoundException e)
			{
				throw new IllegalArgumentException("Failed to open file: " + filename + ".");
			}
			finally
			{
				s.close();
			}
		}

		public String toString() {
			StringBuilder sb = new StringBuilder();
			for(int j = 0; j < rows(); j++) {
				if(j > 0)
					sb.append("\n");
				sb.append(Vec.toString(row(j)));
			}
			return sb.toString();
		}

		/// Saves the matrix to an ARFF file
		public void saveARFF(String filename)
		{		
			PrintWriter os = null;
			
			try
			{
				os = new PrintWriter(filename);
				// Print the relation name, if one has been provided ('x' is default)
				os.print("@RELATION ");
				os.println(m_filename.isEmpty() ? "x" : m_filename);
				
				// Print each attribute in order
				for (int i = 0; i < m_attr_name.size(); i++)
				{
					os.print("@ATTRIBUTE ");
					
					String attributeName = m_attr_name.get(i);
					os.print(attributeName.isEmpty() ? "x" : attributeName);
					
					int vals = valueCount(i);
					
					if (vals == 0) os.println(" REAL");
					else
					{
						os.print(" {");
						for (int j = 0; j < vals; j++)
						{
							os.print(attrValue(i, j));
							if (j + 1 < vals) os.print(",");
						}
						os.println("}");
					}
				}
				
				// Print the data
				os.println("@DATA");
				for (int i = 0; i < rows(); i++)
				{
					double[] row = m_data.get(i);
					for (int j = 0; j < cols(); j++)
					{
						if (row[j] == UNKNOWN_VALUE)
							os.print("?");
						else
						{
							int vals = valueCount(j);
							if (vals == 0) os.print(row[j]);
							else
							{
								int val = (int)row[j];
								if (val >= vals) throw new IllegalArgumentException("Value out of range.");
								os.print(attrValue(j, val));
							}
						}
						
						if (j + 1 < cols())	os.print(",");
					}
					os.println();
				}
			}
			catch (FileNotFoundException e)
			{
				throw new IllegalArgumentException("Error creating file: " + filename + ".");
			}
			finally
			{
				os.close();
			}
		}

		/// Makes a rows-by-columns matrix of *ALL CONTINUOUS VALUES*.
		/// This method wipes out any data currently in the matrix. It also
		/// wipes out any meta-data.
		public void setSize(int rows, int cols)
		{
			m_data.clear();

			// Set the meta-data
			m_filename = "";
			m_attr_name.clear();
			m_str_to_enum.clear();
			m_enum_to_str.clear();

			// Make space for each of the columns, then each of the rows
			newColumns(cols);
			newRows(rows);
		}

		/// Clears this matrix and copies the meta-data from that matrix.
		/// In other words, it makes a zero-row matrix with the same number
		/// of columns as "that" matrix. You will need to call newRow or newRows
		/// to give the matrix some rows.
		@SuppressWarnings("unchecked")
		public void copyMetaData(Matrix that)
		{
			m_data.clear();
			m_attr_name = new ArrayList<String>(that.m_attr_name);
			
			// Make a deep copy of that.m_str_to_enum
			m_str_to_enum = new ArrayList<HashMap<String, Integer>>();
			for (HashMap<String, Integer> map : that.m_str_to_enum)
			{
				HashMap<String, Integer> temp = new HashMap<String, Integer>();
				for (Map.Entry<String, Integer> entry : map.entrySet())
					temp.put(entry.getKey(), entry.getValue());
				
				m_str_to_enum.add(temp);
			}
			
			// Make a deep copy of that.m_enum_to_string
			m_enum_to_str = new ArrayList<HashMap<Integer, String>>();
			for (HashMap<Integer, String> map : that.m_enum_to_str)
			{
				HashMap<Integer, String> temp = new HashMap<Integer, String>();
				for (Map.Entry<Integer, String> entry : map.entrySet())
					temp.put(entry.getKey(), entry.getValue());
				
				m_enum_to_str.add(temp);
			}
		}

		/// Adds a column to this matrix with the specified number of values. (Use 0 for
		/// a continuous attribute.) This method also sets the number of rows to 0, so
		/// you will need to call newRow or newRows when you are done adding columns.
		public void newColumn(int vals)
		{
			m_data.clear();
			String name = "col_" + cols();
			
			m_attr_name.add(name);
			
			HashMap<String, Integer> temp_str_to_enum = new HashMap<String, Integer>();
			HashMap<Integer, String> temp_enum_to_str = new HashMap<Integer, String>();
			
			for (int i = 0; i < vals; i++)
			{
				String sVal = "val_" + i;
				temp_str_to_enum.put(sVal, i);
				temp_enum_to_str.put(i, sVal);
			}
			
			m_str_to_enum.add(temp_str_to_enum);
			m_enum_to_str.add(temp_enum_to_str);
		}
		
		/// Adds a column to this matrix with 0 values (continuous data).
		public void newColumn()
		{
			this.newColumn(0);
		}
		
		/// Adds n columns to this matrix, each with 0 values (continuous data).
		public void newColumns(int n)
		{
			for (int i = 0; i < n; i++)
				newColumn();
		}
		
		/// Adds one new row to this matrix. Returns a reference to the new row.
		public double[] newRow()
		{
			int c = cols();
			if (c == 0)
				throw new IllegalArgumentException("You must add some columns before you add any rows.");
			double[] newRow = new double[c];
			m_data.add(newRow);
			return newRow;
		}
		
		/// Adds 'n' new rows to this matrix
		public void newRows(int n)
		{
			for (int i = 0; i < n; i++)
				newRow();
		}
		
		/// Returns the number of rows in the matrix
		public int rows() { return m_data.size(); }
		
		/// Returns the number of columns (or attributes) in the matrix
		public int cols() { return m_attr_name.size(); }
		
		/// Returns the name of the specified attribute
		public String attrName(int col) { return m_attr_name.get(col); }
		
		/// Returns the name of the specified value
		public String attrValue(int attr, int val)
		{		
			String value = m_enum_to_str.get(attr).get(val);
			if (value == null)
				throw new IllegalArgumentException("No name.");
			else return value;
		}
		
		/// Returns a reference to the specified row
		public double[] row(int index) { return m_data.get(index); }
		
		/// Swaps the positions of the two specified rows
		public void swapRows(int a, int b)
		{
			double[] temp = m_data.get(a);
			m_data.set(a, m_data.get(b));
			m_data.set(b, temp);
		}
		
		/// Returns the number of values associated with the specified attribute (or column)
		/// 0 = continuous, 2 = binary, 3 = trinary, etc.
		public int valueCount(int attr) { return m_enum_to_str.get(attr).size(); }
		
		/// Returns the mean of the elements in the specified column. (Elements with the value UNKNOWN_VALUE are ignored.)
		public double columnMean(int col)
		{
			double sum = 0.0;
			int count = 0;
			for (double[] list : m_data)
			{
				double val = list[col];
				if (val != UNKNOWN_VALUE)
				{
					sum += val;
					count++;
				}
			}
			
			return sum / count;
		}
		
		/// Returns the minimum element in the specified column. (Elements with the value UNKNOWN_VALUE are ignored.)
		public double columnMin(int col)
		{
			double min = Double.MAX_VALUE;
			for (double[] list : m_data)
			{
				double val = list[col];
				if (val != UNKNOWN_VALUE)
					min = Math.min(min, val);
			}
			
			return min;
		}

		/// Returns the maximum element in the specifed column. (Elements with the value UNKNOWN_VALUE are ignored.)
		public double columnMax(int col)
		{
			double max = -Double.MAX_VALUE;
			for (double[] list : m_data)
			{
				double val = list[col];
				if (val != UNKNOWN_VALUE)
					max = Math.max(max, val);
			}
			
			return max;
		}
		
		/// Returns the most common value in the specified column. (Elements with the value UNKNOWN_VALUE are ignored.)
		public double mostCommonValue(int col)
		{
			HashMap<Double, Integer> counts = new HashMap<Double, Integer>();
			for (double[] list : m_data)
			{
				double val = list[col];
				if (val != UNKNOWN_VALUE)
				{
					Integer result = counts.get(val);
					if (result == null) result = 0;
					
					counts.put(val, result + 1);
				}
			}
			
			int valueCount = 0;
			double value   = 0;
			for (Map.Entry<Double, Integer> entry : counts.entrySet())
			{
				if (entry.getValue() > valueCount)
				{
					value      = entry.getKey();
					valueCount = entry.getValue();
				}
			}
			
			return value;
		}

		/// Copies the specified rectangular portion of that matrix, and puts it in the specified location in this matrix.
		public void copyBlock(int destRow, int destCol, Matrix that, int rowBegin, int colBegin, int rowCount, int colCount)
		{
			if (destRow + rowCount > this.rows() || destCol + colCount > this.cols())
				throw new IllegalArgumentException("Out of range for destination matrix.");
			if (rowBegin + rowCount > that.rows() || colBegin + colCount > that.cols())
				throw new IllegalArgumentException("Out of range for source matrix.");

			// Copy the specified region of meta-data
			for (int i = 0; i < colCount; i++)
			{
				m_attr_name.set(destCol + i, that.m_attr_name.get(colBegin + i));
				m_str_to_enum.set(destCol + i, new HashMap<String, Integer>(that.m_str_to_enum.get(colBegin + i)));
				m_enum_to_str.set(destCol + i, new HashMap<Integer, String>(that.m_enum_to_str.get(colBegin + i)));
			}

			// Copy the specified region of data
			for (int i = 0; i < rowCount; i++)
			{
				double[] source = that.row(rowBegin + i);
				double[] dest = this.row(destRow + i);
				for(int j = 0; j < colCount; j++)
					dest[j] = source[colBegin + j];
			}
		}
		
		/// Sets every element in the matrix to the specified value.
		public void setAll(double val)
		{
			for (double[] list : m_data) {
				for(int i = 0; i < list.length; i++)
					list[i] = val;
			}
		}

		/// Sets this to the identity matrix.
		public void setToIdentity()
		{
			setAll(0.0);
			int m = Math.min(cols(), rows());
			for(int i = 0; i < m; i++)
				row(i)[i] = 1.0;
		}

		/// Throws an exception if that has a different number of columns than
		/// this, or if one of its columns has a different number of values.
		public void checkCompatibility(Matrix that)
		{
			int c = cols();
			if (that.cols() != c)
				throw new IllegalArgumentException("Matrices have different number of columns.");
			
			for (int i = 0; i < c; i++)
			{
				if (valueCount(i) != that.valueCount(i))
					throw new IllegalArgumentException("Column " + i + " has mis-matching number of values.");
			}
		}
	}


	public static class LayerTanh {
		public Matrix weights; // rows are inputs, cols are outputs
		public double[] bias;
		public double[] net;
		public double[] activation;
		public double[] error;


		LayerTanh(int inputs, int outputs) {
			weights = new Matrix();
			weights.setSize(inputs, outputs);
			bias = new double[outputs];
			net = new double[outputs];
			activation = new double[outputs];
			error = new double[outputs];
		}


		LayerTanh(LayerTanh that) {
			weights = new Matrix(that.weights);
			bias = Vec.copy(that.bias);
			net = Vec.copy(that.net);
			activation = Vec.copy(that.activation);
			error = Vec.copy(that.error);
		}


		void copy(LayerTanh src) {
			if(src.weights.rows() != weights.rows() || src.weights.cols() != weights.cols())
				throw new IllegalArgumentException("mismatching sizes");
			weights.setSize(src.weights.rows(), src.weights.cols());
			weights.copyBlock(0, 0, src.weights, 0, 0, src.weights.rows(), src.weights.cols());
			for(int i = 0; i < bias.length; i++) {
				bias[i] = src.bias[i];
			}
		}


		int inputCount() { return weights.rows(); }
		int outputCount() { return weights.cols(); }


		void initWeights(Random r) {
			double dev = Math.max(0.3, 1.0 / weights.rows());
			for(int i = 0; i < weights.rows(); i++) {
				double[] row = weights.row(i);
				for(int j = 0; j < weights.cols(); j++) {
					row[j] = dev * r.nextGaussian();
				}
			}
			for(int j = 0; j < weights.cols(); j++) {
				bias[j] = dev * r.nextGaussian();
			}
		}


		int countWeights() {
			return weights.rows() * weights.cols() + bias.length;
		}


		int setWeights(double[] w, int start) {
			int oldStart = start;
			for(int i = 0; i < bias.length; i++)
				bias[i] = w[start++];
			for(int i = 0; i < weights.rows(); i++)
			{
				double[] row = weights.row(i);
				for(int j = 0; j < weights.cols(); j++)
					row[j] = w[start++];
			}
			return start - oldStart;
		}


		void feedForward(double[] in) {
			if(in.length != weights.rows())
				throw new IllegalArgumentException("size mismatch. " + Integer.toString(in.length) + " != " + Integer.toString(weights.rows()));
			for(int i = 0; i < net.length; i++)
				net[i] = bias[i];
			for(int j = 0; j < weights.rows(); j++) {
				double v = in[j];
				double[] w = weights.row(j);
				for(int i = 0; i < weights.cols(); i++)
					net[i] += v * w[i];
			}
		}


		void feedForward2(double[] in1, double[] in2) {
			if(in1.length + in2.length != weights.rows())
				throw new IllegalArgumentException("size mismatch. " + Integer.toString(in1.length) + " + " + Integer.toString(in2.length) + " != " + Integer.toString(weights.rows()));
			for(int i = 0; i < net.length; i++)
				net[i] = bias[i];
			for(int j = 0; j < in1.length; j++) {
				double v = in1[j];
				double[] w = weights.row(j);
				for(int i = 0; i < weights.cols(); i++)
					net[i] += v * w[i];
			}
			for(int j = 0; j < in2.length; j++) {
				double v = in2[j];
				double[] w = weights.row(in1.length + j);
				for(int i = 0; i < weights.cols(); i++)
					net[i] += v * w[i];
			}
		}


		void activate() {
			for(int i = 0; i < net.length; i++) {
				activation[i] = Math.tanh(net[i]);
			}
		}


		void computeError(double[] target) {
			if(target.length != activation.length)
				throw new IllegalArgumentException("size mismatch. " + Integer.toString(target.length) + " != " + Integer.toString(activation.length));
			for(int i = 0; i < activation.length; i++) {
				if(target[i] < -1.0 || target[i] > 1.0)
					throw new IllegalArgumentException("target value out of range for the tanh activation function");
				error[i] = target[i] - activation[i];
			}
		}


		void deactivate() {
			for(int i = 0; i < error.length; i++) {
				error[i] *= (1.0 - activation[i] * activation[i]);
			}
		}


		void feedBack(double[] upstream) {
			if(upstream.length != weights.rows())
				throw new IllegalArgumentException("size mismatch");
			for(int j = 0; j < weights.rows(); j++) {
				double[] w = weights.row(j);
				double d = 0.0;
				for(int i = 0; i < weights.cols(); i++) {
					d += error[i] * w[i];
				}
				upstream[j] = d;
			}
		}


		void refineInputs(double[] inputs, double learningRate) {
			if(inputs.length != weights.rows())
				throw new IllegalArgumentException("size mismatch");
			for(int j = 0; j < weights.rows(); j++) {
				double[] w = weights.row(j);
				double d = 0.0;
				for(int i = 0; i < weights.cols(); i++) {
					d += error[i] * w[i];
				}
				inputs[j] += learningRate * d;
			}
		}


		void updateWeights(double[] in, double learningRate) {
			for(int i = 0; i < bias.length; i++) {
				bias[i] += learningRate * error[i];
			}
			for(int j = 0; j < weights.rows(); j++) {
				double[] w = weights.row(j);
				double x = learningRate * in[j];
				for(int i = 0; i < weights.cols(); i++) {
					w[i] += x * error[i];
				}
			}
		}

		// Applies both L2 and L1 regularization to the weights and bias values
		void regularizeWeights(double lambda) {
			for(int i = 0; i < weights.rows(); i++) {
				double[] row = weights.row(i);
				for(int j = 0; j < row.length; j++) {
					row[j] *= (1.0 - lambda);
					if(row[j] < 0.0)
						row[j] += lambda;
					else
						row[j] -= lambda;
				}
			}
			for(int j = 0; j < bias.length; j++) {
				bias[j] *= (1.0 - lambda);
				if(bias[j] < 0.0)
					bias[j] += lambda;
				else
					bias[j] -= lambda;
			}
		}
	}




	public static class NeuralNet {
		public ArrayList<LayerTanh> layers;


		/// General-purpose constructor. (Starts with no layers. You must add at least one.)
		NeuralNet() {
			layers = new ArrayList<LayerTanh>();
		}


		/// Copy constructor
		NeuralNet(NeuralNet that) {
			layers = new ArrayList<LayerTanh>();
			for(int i = 0; i < that.layers.size(); i++) {
				layers.add(new LayerTanh(that.layers.get(i)));
			}
		}


		/// Initializes the weights and biases with small random values
		void init(Random r) {
			for(int i = 0; i < layers.size(); i++) {
				layers.get(i).initWeights(r);
			}
		}


		/// Copies all the weights and biases from "that" into "this".
		/// (Assumes the corresponding topologies already match.)
		void copy(NeuralNet that) {
			if(layers.size() != that.layers.size())
				throw new IllegalArgumentException("Unexpected number of layers");
			for(int i = 0; i < layers.size(); i++) {
				layers.get(i).copy(that.layers.get(i));
			}
		}


		/// Feeds "in" into this neural network and propagates it forward to compute predicted outputs.
		double[] forwardProp(double[] in) {
			LayerTanh l = null;
			for(int i = 0; i < layers.size(); i++) {
				l = layers.get(i);
				l.feedForward(in);
				l.activate();
				in = l.activation;
			}
			return l.activation;
		}


		/// Feeds the concatenation of "in1" and "in2" into this neural network and propagates it forward to compute predicted outputs.
		double[] forwardProp2(double[] in1, double[] in2) {
			LayerTanh l = layers.get(0);
			l.feedForward2(in1, in2);
			l.activate();
			double[] in = l.activation;
			for(int i = 1; i < layers.size(); i++) {
				l = layers.get(i);
				l.feedForward(in);
				l.activate();
				in = l.activation;
			}
			return l.activation;
		}


		/// Backpropagates the error to the upstream layer.
		void backProp(double[] target) {
			int i = layers.size() - 1;
			LayerTanh l = layers.get(i);
			l.computeError(target);
			l.deactivate();
			for(i--; i >= 0; i--) {
				LayerTanh upstream = layers.get(i);
				l.feedBack(upstream.error);
				upstream.deactivate();
				l = upstream;
			}
		}


		/// Backpropagates the error from another neural network. (This is used when training autoencoders.)
		void backPropFromDecoder(NeuralNet decoder) {
			int i = layers.size() - 1;
			LayerTanh l = decoder.layers.get(0);
			LayerTanh upstream = layers.get(i);
			l.feedBack(upstream.error);
			l = upstream;
			//l.bendHinge(learningRate);
			l.deactivate();
			for(i--; i >= 0; i--) {
				upstream = layers.get(i);
				l.feedBack(upstream.error);
				//upstream.bendHinge(learningRate);
				upstream.deactivate();
				l = upstream;
			}
		}


		/// Updates the weights and biases
		void descendGradient(double[] in, double learningRate) {
			for(int i = 0; i < layers.size(); i++) {
				LayerTanh l = layers.get(i);
				l.updateWeights(in, learningRate);
				in = l.activation;
			}
		}


		/// Keeps the weights and biases from getting too big
		void regularize(double learningRate, double lambda) {
			double amount = learningRate * lambda;
			double smallerAmount = 0.1 * amount;
			for(int i = 0; i < layers.size(); i++) {
				LayerTanh lay = layers.get(i);
				//lay.straightenHinge(amount);
				lay.regularizeWeights(smallerAmount);
			}
		}


		/// Refines the weights and biases with on iteration of stochastic gradient descent.
		void trainIncremental(double[] in, double[] target, double learningRate) {
			forwardProp(in);
			backProp(target);
			//backPropAndBendHinge(target, learningRate);
			descendGradient(in, learningRate);
		}


		/// Refines "in" with one iteration of stochastic gradient descent.
		void refineInputs(double[] in, double[] target, double learningRate) {
			forwardProp(in);
			backProp(target);
			layers.get(0).refineInputs(in, learningRate);
		}


		static void testMath() {
			NeuralNet nn = new NeuralNet();
			LayerTanh l1 = new LayerTanh(2, 3);
			l1.weights.row(0)[0] = 0.1;
			l1.weights.row(0)[1] = 0.0;
			l1.weights.row(0)[2] = 0.1;
			l1.weights.row(1)[0] = 0.1;
			l1.weights.row(1)[1] = 0.0;
			l1.weights.row(1)[2] = -0.1;
			l1.bias[0] = 0.1;
			l1.bias[1] = 0.1;
			l1.bias[2] = 0.0;
			nn.layers.add(l1);

			LayerTanh l2 = new LayerTanh(3, 2);
			l2.weights.row(0)[0] = 0.1;
			l2.weights.row(0)[1] = 0.1;
			l2.weights.row(1)[0] = 0.1;
			l2.weights.row(1)[1] = 0.3;
			l2.weights.row(2)[0] = 0.1;
			l2.weights.row(2)[1] = -0.1;
			l2.bias[0] = 0.1;
			l2.bias[1] = -0.2;
			nn.layers.add(l2);

			System.out.println("l1 weights:" + l1.weights.toString());
			System.out.println("l1 bias:" + Vec.toString(l1.bias));
			System.out.println("l2 weights:" + l2.weights.toString());
			System.out.println("l2 bias:" + Vec.toString(l2.bias));

			System.out.println("----Forward prop");
			double in[] = new double[2];
			in[0] = 0.3;
			in[1] = -0.2;
			double[] out = nn.forwardProp(in);
			System.out.println("activation:" + Vec.toString(out));

			System.out.println("----Back prop");
			double targ[] = new double[2];
			targ[0] = 0.1;
			targ[1] = 0.0;
			nn.backProp(targ);
			System.out.println("error 2:" + Vec.toString(l2.error));
			System.out.println("error 1:" + Vec.toString(l1.error));
			
			nn.descendGradient(in, 0.1);
			System.out.println("----Descending gradient");
			System.out.println("l1 weights:" + l1.weights.toString());
			System.out.println("l1 bias:" + Vec.toString(l1.bias));
			System.out.println("l2 weights:" + l2.weights.toString());
			System.out.println("l2 bias:" + Vec.toString(l2.bias));

			if(Math.abs(l1.weights.row(0)[0] - 0.10039573704287) > 0.0000000001)
				throw new IllegalArgumentException("failed");
			if(Math.abs(l1.weights.row(0)[1] - 0.0013373814241446) > 0.0000000001)
				throw new IllegalArgumentException("failed");
			if(Math.abs(l1.bias[1] - 0.10445793808048) > 0.0000000001)
				throw new IllegalArgumentException("failed");
			System.out.println("passed");
		}

		public static void testVisual() throws Exception {
			// Make some data
			Random rand = new Random(1234);
			Matrix features = new Matrix();
			features.setSize(1000, 2);
			Matrix labels = new Matrix();
			labels.setSize(1000, 2);
			for(int i = 0; i < 1000; i++) {
				
				double x = rand.nextDouble() * 2 - 1;
				double y = rand.nextDouble() * 2 - 1;
				features.row(i)[0] = x;
				features.row(i)[1] = y;
				labels.row(i)[0] = (y < x * x ? 0.9 : 0.1);
				labels.row(i)[1] = (x < y * y ? 0.1 : 0.9);
			}

			// Train on it
			NeuralNet nn = new NeuralNet();
			nn.layers.add(new LayerTanh(2, 30));
			nn.layers.add(new LayerTanh(30, 2));
			nn.init(rand);
			int iters = 10000000;
			double learningRate = 0.01;
			double lambda = 0.0001;
			for(int i = 0; i < iters; i++) {
				int index = rand.nextInt(features.rows());
				nn.regularize(learningRate, lambda);
				nn.trainIncremental(features.row(index), labels.row(index), 0.01);
				if(i % 1000000 == 0)
					System.out.println(Double.toString(((double)i * 100)/ iters) + "%");
			}

			// Visualize it
			for(int i = 0; i < nn.layers.size(); i++) {
				System.out.print("Layer " + Integer.toString(i) + ": ");
//				for(int j = 0; j < nn.layers.get(i).hinge.length; j++)
//					System.out.print(Double.toString(nn.layers.get(i).hinge[j]) + ", ");
				System.out.println();
			}
			BufferedImage image = new BufferedImage(100, 200, BufferedImage.TYPE_INT_ARGB);
			double[] in = new double[2];
			for(int y = 0; y < 100; y++) {
				for(int x = 0; x < 100; x++) {
					in[0] = ((double)x) / 100 * 2 - 1;
					in[1] = ((double)y) / 100 * 2 - 1;
					double[] out = nn.forwardProp(in);
					int g = Math.max(0, Math.min(255, (int)(out[0] * 256)));
					image.setRGB(x, y, new Color(g, g, g).getRGB());
					g = Math.max(0, Math.min(255, (int)(out[1] * 256)));
					image.setRGB(x, y + 100, new Color(g, g, g).getRGB());
				}
			}
			ImageIO.write(image, "png", new File("viz.png"));
		}
	}

	
	int index; // a temporary value used to pass values around
	NeuralNet nn;
	double[] in;

	BellBrent2() {

	}

	BellBrent2(double[] weights) {
		in = new double[20];
		nn = new NeuralNet();
		nn.layers.add(new LayerTanh(in.length, 8));
		nn.layers.add(new LayerTanh(8, 10));
		nn.layers.add(new LayerTanh(10, 3));
		setWeights(weights);
	}
	public void reset() {
		double[] finalGame = {0.21112732349160243, -1.3284286154538962, -0.5477090868688381, 0.8362676057663732, 1.6178594410191631, 0.09883730060853195, 0.10612072232413444, -0.4952776967986526, 0.29762938907675124, -0.3829582880859058, 0.48441133804688674, 0.22935434191744627, 0.5507866634091275, 0.6514415422603119, 0.9168547666057252, 0.08575736642975756, -1.4029022265050661, 0.19304220377143255, -0.023626603076140066, -2.145144715773424, 1.2551006540262637, -0.38370871859309613, 0.8963331096474713, -0.4048098780745281, 1.197403424546195, -0.907827879369115, 2.2732040633549992, 0.02077670321363306, 1.0244578266076074, -0.42056670008849606, 0.45543866482101375, -2.0311815298932108, -1.625613802990942, -0.8842981833441458, 1.2480321927945583, -0.15952617075720013, -1.1431351484266457, -1.2704313130988962, -0.49530974609634515, 0.1620353058631801, -0.15848534622130522, -0.6207006925170612, -0.115439738457424, -1.4591108541724098, 0.46559973106718505, -0.8066199153257931, -0.12267354175955487, 0.3065223023697042, -0.15784855677859544, -1.3921227541853445, -0.7574309879284145, 0.5165480927339023, 0.41219367887355757, -0.11508290918914793, -0.9840893304299225, -2.345129044699004, -0.17162842931121167, -1.352795648955529, 0.14250280051822437, 0.45375031978694896, -0.06947134604837545, 0.3616841209742882, 1.8692538615170224, 0.03281612527338388, -0.7981391645227384, 0.1379640244813684, -1.2452054255917553, -0.5294719501984457, -0.8007103735554741, 0.1553367055556964, -0.30138411312595975, -0.19658901368127296, 1.4213075290785573, 0.3298542202849766, -0.9454066422388864, -0.5077788439460371, -0.1649593028641593, 0.03859015014332201, -0.06338992911867596, 0.9336412952298779, -1.3592892946460131, -0.21841348552167625, -0.17008526136281918, 0.17495392853500824, 0.44802909316271644, 0.46701045943726266, -0.016328781933428932, -1.1018558564857313, 0.2126829739474869, -1.005855889955682, 1.142263296842644, 0.89027212739276, -1.6185625149157228, 0.31524068784500636, 1.2024774211034885, -0.14936063540186892, 0.11939229737467817, -1.806515612606433, -0.5794211722708863, -0.20215021369862934, -0.6448404025098142, -0.6130188643687791, 0.19064962502346824, -0.5477670722732061, -1.1522466333805201, -1.4087976439493983, -0.13506537720854128, 1.2475873319687658, 0.17979980535716167, -1.9846147545865835, 0.299507255454186, -0.6980381512194145, -2.2097940960778053, 0.025324607437965684, 1.6596437457781759, -1.3186315590370072, 0.1762574794161522, -1.5551532073324423, 0.6676930601875634, 0.9454816256556868, -0.23501841256421596, -1.4233951539782472, 0.4897232804519044, -1.8091657156730196, 1.763024714488278, 3.0874756115365653E-4, -2.0811858303610045, -0.3918366860225748, -1.3251802572502163, -0.5459969154012323, 0.32301754912240693, -1.032135992124288, 0.6284684139892894, 0.36531748605339787, -1.35859960778247, -1.0748408793048014, -2.1728918517422597, 1.001222665750878, 0.35566385049233973, -1.5156395358111725, 2.1268235207548174, -0.4370672779112322, 0.8925213653707742, -0.45759679247313995, -0.49349650000925377, 0.36427652802644594, -0.6588271995772423, 0.2605559641583875, -0.22883141131691886, 0.2539547754220576, -1.5711503097181188, 1.3522393946595503, -0.9356222612114199, -0.9855619645024135, -1.1451558734667617, 0.7807244165801822, -0.22395266324319263, 0.8854687919573919, 0.7930834373427578, -1.05189493087308, -1.5007639757915567, -0.22874198156794434, 1.3739727961916914, 0.8266862762180158, 0.0902414090449021, -1.3515702787586414, 0.39684789283805116, -1.1854216770999444, 1.5509284272908983, 1.7499562989611543, 0.42469847147060813, -0.34026923070332327, -2.798885531700293, -0.45928126950763437, -1.3452496163270622, 3.2242710646992205, 0.022458500192526463, -0.8789033569840394, 0.09708379795283577, -0.5459443388171975, -0.4579471792594911, -1.0289195467748875, 0.17822140692997895, 1.6114257456788166, -0.038369406469679734, -1.2442857049895073, -0.1338688004903728, -0.13897137034103327, 0.7106738888699428, -1.4150824055668356, 1.5363382123692142, 0.34861286568115646, -0.4672110997349924, 0.5321660436396598, -2.1150199222688797, -1.7893427646112479, 0.74712758232522, 2.463027928012655, -0.6569237117670063, 0.6521988505503274, 1.57141675671424, 0.5145983187379328, 0.6675401912382697, 0.0382558454734811, -1.376417557164395, -1.7853292558993588, -0.03537756937153658, 0.3485201830347864, -0.6433733126296574, -0.381981387782906, 0.9884667404524718, 1.7960770402496586, 0.6293388197777837, -0.9307977239777632, 1.7390443071570245, 0.4203242648270653, -0.21974087843993592, 1.6882833037760479, -0.2985147904656494, -1.008752375588371, -1.7707629813199297, -0.17073120361473323, -0.1943328577635861, -0.3992599376254279, 1.735505921551561, 0.1082493954861698, 0.7817167991033996, 0.1940729420365212, 1.993232627349531, -1.6696397659746642, 0.2597054783551031, -1.629832142528843, 0.1301564558845534, 0.6026882760454452, -0.45883360576795207, 0.13247100536874357, -1.5466837050853661, 0.5299090966232399, 0.5044182136790847, -1.366770578591991, 1.8751097702999344, 0.6607005472956359, -0.24307920333245678, 1.4464091466788047, 0.39980525241032105, 3.3142298219453608, -0.05282977757289053, 0.4092534864536279, 1.7420120466735614, -1.1651625164111068, -0.1299365952376272, -1.2310741018748002, 1.2195494944996301, -0.3968505033217308, -1.1754544760350925, 0.47202333596098534, 1.784243461333044, -1.5795735226963823, 1.215007970997528, -0.9266452128139055, 1.0687385872012922, -2.888979164053574, 0.8034527336036166, -1.3511946348775992, -1.897573373423779, 0.8143261409694597, 0.5286261347849008, -0.44718894487518585, 0.17015679590418983, -0.6472075596560066, -2.259247612936252, -0.95695047251063, -0.4059503563959095, -1.431194657887998, -0.2166407085852435, -0.4303915891482964, -0.8832471381651369, -2.018215775963517, 0.7631648072039497, -1.0650401056104295, -1.2288893184463108, -1.0943492468572953, -1.9108653674113525, -1.0759328152572207, -1.5025575630354893, -1.4891420871332803, -0.47722295663896613, -3.47994576067425, 0.6274393912015482, 1.7879457883897458, 0.30699198141424927};
		in = new double[20];
		nn = new NeuralNet();
		nn.layers.add(new LayerTanh(in.length, 8));
		nn.layers.add(new LayerTanh(8, 10));
		nn.layers.add(new LayerTanh(10, 3));
		setWeights(finalGame);
	}

	/// Returns the number of weights necessary to fully-parameterize this agent
	int countWeights() {
		int n = 0;
		for(int i = 0; i < nn.layers.size(); i++)
			n += nn.layers.get(i).countWeights();
		return n;
	}


	/// Sets the parameters of this agent with the specified weights
	void setWeights(double[] weights) {
		if(weights.length != countWeights())
			throw new IllegalArgumentException("Wrong number of weights. Got " + Integer.toString(weights.length) + ", expected " + Integer.toString(countWeights()));
		int start = 0;
		for(int i = 0; i < nn.layers.size(); i++)
			start += nn.layers.get(i).setWeights(weights, start);
	}


	public static float sq_dist(float x1, float y1, float x2, float y2) {
		return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
	}

	float nearestBombTarget(Model m, float x, float y) {
		index = -1;
		float dd = Float.MAX_VALUE;
		for(int i = 0; i < m.getBombCount(); i++) {
			float d = sq_dist(x, y, m.getBombTargetX(i), m.getBombTargetY(i));
			if(d < dd) {
				dd = d;
				index = i;
			}
		}
		return dd;
	}

	float nearestOpponent(Model m, float x, float y) {
		index = -1;
		float dd = Float.MAX_VALUE;
		for(int i = 0; i < m.getSpriteCountOpponent(); i++) {
			if(m.getEnergyOpponent(i) < 0)
				continue; // don't care about dead opponents
			float d = sq_dist(x, y, m.getXOpponent(i), m.getYOpponent(i));
			if(d < dd) {
				dd = d;
				index = i;
			}
		}
		return dd;
	}

	void avoidBombs(Model m, int i) {
		if(nearestBombTarget(m, m.getX(i), m.getY(i)) <= 2.0f * Model.BLAST_RADIUS * Model.BLAST_RADIUS) {
			float dx = m.getX(i) - m.getBombTargetX(index);
			float dy = m.getY(i) - m.getBombTargetY(index);
			if(dx == 0 && dy == 0)
				dx = 1.0f;
			m.setDestination(i, m.getX(i) + dx * 10.0f, m.getY(i) + dy * 10.0f);
		}
	}

	void beDefender(Model m, int i) {
		// Find the opponent nearest to my flag
		nearestOpponent(m, Model.XFLAG, Model.YFLAG);
		if(index >= 0) {
			float enemyX = m.getXOpponent(index);
			float enemyY = m.getYOpponent(index);

			// Stay between the enemy and my flag
			m.setDestination(i, 0.5f * (Model.XFLAG + enemyX), 0.5f * (Model.YFLAG + enemyY));

			// Throw boms if the enemy gets close enough
			if(sq_dist(enemyX, enemyY, m.getX(i), m.getY(i)) <= Model.MAX_THROW_RADIUS * Model.MAX_THROW_RADIUS)
				m.throwBomb(i, enemyX, enemyY);
		}
		else {
			// Guard the flag
			m.setDestination(i, Model.XFLAG + Model.MAX_THROW_RADIUS, Model.YFLAG);
		}

		// If I don't have enough energy to throw a bomb, rest
		if(m.getEnergySelf(i) < Model.BOMB_COST)
			m.setDestination(i, m.getX(i), m.getY(i));

		// Try not to die
		avoidBombs(m, i);
	}

	void beFlagAttacker(Model m, int i) {
		// Head for the opponent's flag
		m.setDestination(i, Model.XFLAG_OPPONENT - Model.MAX_THROW_RADIUS + 1, Model.YFLAG_OPPONENT);

		// Shoot at the flag if I can hit it
		if(sq_dist(m.getX(i), m.getY(i), Model.XFLAG_OPPONENT, Model.YFLAG_OPPONENT) <= Model.MAX_THROW_RADIUS * Model.MAX_THROW_RADIUS) {
			m.throwBomb(i, Model.XFLAG_OPPONENT, Model.YFLAG_OPPONENT);
		}

		// Try not to die
		avoidBombs(m, i);
	}

	void beAggressor(Model m, int i) {
		float myX = m.getX(i);
		float myY = m.getY(i);

		// Find the opponent nearest to me
		nearestOpponent(m, myX, myY);
		if(index >= 0) {
			float enemyX = m.getXOpponent(index);
			float enemyY = m.getYOpponent(index);

			if(m.getEnergySelf(i) >= m.getEnergyOpponent(index)) {

				// Get close enough to throw a bomb at the enemy
				float dx = myX - enemyX;
				float dy = myY - enemyY;
				float t = 1.0f / Math.max(Model.EPSILON, (float)Math.sqrt(dx * dx + dy * dy));
				dx *= t;
				dy *= t;
				m.setDestination(i, enemyX + dx * (Model.MAX_THROW_RADIUS - Model.EPSILON), enemyY + dy * (Model.MAX_THROW_RADIUS - Model.EPSILON));

				// Throw bombs
				if(sq_dist(enemyX, enemyY, m.getX(i), m.getY(i)) <= Model.MAX_THROW_RADIUS * Model.MAX_THROW_RADIUS)
					m.throwBomb(i, enemyX, enemyY);
			}
			else {

				// If the opponent is close enough to shoot at me...
				if(sq_dist(enemyX, enemyY, myX, myY) <= (Model.MAX_THROW_RADIUS + Model.BLAST_RADIUS) * (Model.MAX_THROW_RADIUS + Model.BLAST_RADIUS)) {
					m.setDestination(i, myX + 10.0f * (myX - enemyX), myY + 10.0f * (myY - enemyY)); // Flee
				}
				else {
					m.setDestination(i, myX, myY); // Rest
				}
			}
		}

		// Try not to die
		avoidBombs(m, i);
	}

	public void update(Model m) {

		// Compute some features
		in[0] = m.getX(0) / 600.0 - 0.5;
		in[1] = m.getY(0) / 600.0 - 0.5;
		in[2] = m.getX(1) / 600.0 - 0.5;
		in[3] = m.getY(1) / 600.0 - 0.5;
		in[4] = m.getX(2) / 600.0 - 0.5;
		in[5] = m.getY(2) / 600.0 - 0.5;
		in[6] = nearestOpponent(m, m.getX(0), m.getY(0)) / 600.0 - 0.5;
		in[7] = nearestOpponent(m, m.getX(0), m.getY(0)) / 600.0 - 0.5;
		in[8] = nearestOpponent(m, m.getX(0), m.getY(0)) / 600.0 - 0.5;
		in[9] = nearestBombTarget(m, m.getX(0), m.getY(0)) / 600.0 - 0.5;
		in[10] = nearestBombTarget(m, m.getX(0), m.getY(0)) / 600.0 - 0.5;
		in[11] = nearestBombTarget(m, m.getX(0), m.getY(0)) / 600.0 - 0.5;
		in[12] = m.getEnergySelf(0);
		in[13] = m.getEnergySelf(1);
		in[14] = m.getEnergySelf(2);
		in[15] = m.getEnergyOpponent(0);
		in[16] = m.getEnergyOpponent(1);
		in[17] = m.getEnergyOpponent(2);
		in[18] = m.getFlagEnergySelf();
		in[19] = m.getFlagEnergyOpponent();

		// Determine what each agent should do
		double[] out = nn.forwardProp(in);

		// Do it
		for(int i = 0; i < 3; i++)
		{
			if(out[i] < -0.333)
				beDefender(m, i);
			else if(out[i] > 0.333)
				beAggressor(m, i);
			else
				beFlagAttacker(m, i);
		}
	}
}
