import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.Iterator;
import java.lang.StringBuilder;

class Game {
	
	private static void Mutate(Matrix population)
	{
		Random r = new Random();
		for(int i = 0; i < 100; i++)
		{
			if(r.nextDouble() < population.row(i)[291])
			{
				int index = Math.abs(r.nextInt() % population.row(i).length);
				//System.out.println("Before" + population.row(i)[index]);
				population.row(i)[index] += r.nextGaussian() * population.row(i)[292];
				//System.out.println("After" + population.row(i)[index]);
				//System.out.println("stddev" + population.row(i)[292]);
				

			}
		}
		for(int i = 0; i < 100; i++)
		{
			
			if(population.row(i)[291] > 0.9)
				population.row(i)[291] = 0.9;
			if(population.row(i)[291] < 0.1)
				population.row(i)[291] = 0.1;
			if(population.row(i)[292] > 10)
				population.row(i)[292] = 10;
			if(population.row(i)[292] < 0.1)
				population.row(i)[292] = 0.1;
			if(population.row(i)[293] > 0.9)
				population.row(i)[293] = 0.9;
			if(population.row(i)[293] < 0.1)
				population.row(i)[293] = 0.1;
			if(population.row(i)[294] > 10)
				population.row(i)[294] = 10;
			if(population.row(i)[294] < 1)
				population.row(i)[294] = 1;
			if(population.row(i)[295] > 20)
				population.row(i)[295] = 20;
			if(population.row(i)[295] < 1)
				population.row(i)[295] = 1;
			if(population.row(i)[296] > 5)
				population.row(i)[296] = 5;
			if(population.row(i)[296] < 1)
				population.row(i)[296] = 1;
		}
	}
	
	private static int[] NaturalSelection(Matrix population) throws Exception
	{
		int[] dead = new int[(int)population.row(0)[294]];
		Random r = new Random();
		for(int i = 0; i < dead.length; i++)
		{
			int agent1Choice = Math.abs(r.nextInt()) % population.rows();
			int agent2Choice = Math.abs(r.nextInt()) % population.rows();
			int winner = 0;
			//System.out.println(population.rows());
			double[] agent1weights = new double[291];
			double[] agent2weights = new double[291];
			for(int temp = 0; temp < 291; temp++)
			{
				agent1weights[temp] = population.row(agent1Choice)[temp];
				agent2weights[temp] = population.row(agent2Choice)[temp];
			}
			//NeuralAgent agent1 = new NeuralAgent(agent1weights);
			//NeuralAgent agent2 = new NeuralAgent(agent2weights);
			try{
				int chooseOpponent = r.nextInt(4);
				IAgent agent = new Winner2016a();
				if(chooseOpponent == 0)
					agent = new Winner2016a();
				else if(chooseOpponent == 1)
					agent = new Winner2015a();
				else if(chooseOpponent == 2)
					agent = new Winner2015b();
				else
					agent = new BellBrent();
				if(r.nextDouble() < 0.5)
					winner = Controller.doBattleNoGui(new BellBrent2(agent1weights), agent);
				else
					winner = Controller.doBattleNoGui(agent, new BellBrent2(agent1weights));					
			}catch(Exception e)
			{
				System.out.println("Exception Thrown:\n" + e.getMessage());
			}
			
			double survive = r.nextDouble();
			if(winner == 1)
			{
				if(survive < population.row(i)[293])
					dead[i] = agent2Choice;
				else
					dead[i] = agent1Choice;	
			}
			else if(winner == -1)
			{
				if(survive < population.row(i)[293])
					dead[i] = agent1Choice;
				else
					dead[i] = agent2Choice;
			}
				
		}
		return dead;
	}
	
	public static void Replenish(Matrix population, int[] dead)
	{
		for(int i = 0; i < dead.length; i++)
		{
			Random r = new Random();
			//Find parent1 and several parent2
			boolean isSame = true;
			int parent1 = Math.abs(r.nextInt()) % population.rows();
			while(isSame)
			{
				parent1 = Math.abs(r.nextInt()) % population.rows();
				isSame = false;
				for(int j = 0; j < dead.length; j++)
				{
					if(dead[j] == parent1)
						isSame = true;
				}
			}
			int[] parent2 = new int[(int)population.row(i)[295]];
			for(int j = 0; j < parent2.length; j++)
			{	
				boolean anySame = true;
				while(anySame) //Keepgoing until all are unique
				{
					int newIndex =  Math.abs(r.nextInt()) % population.rows();
					anySame = false;
					for(int k = 0; k < j; k++)
					{
						for(int l = 0; l < dead.length; l++)
						{
							if(newIndex == parent2[k] || newIndex == dead[l])
								anySame = true;
						}
					}
					parent2[j] = newIndex;
				}
				//System.out.println("Indexes: " + parent2[j]);
			}
			//System.out.println("");
			
			
			//Find Best Match
			double bestMatch = Double.MAX_VALUE;
			int bestIndex = 0;
			for(int j = 0; j < parent2.length; j++)
			{
				double betterMatch = 0;
				for(int k = 0; k < population.cols(); k++)
				{
					betterMatch += Math.pow(Math.abs(Math.pow(population.row(parent1)[k], population.row(parent1)[296]) - Math.pow(population.row(parent2[j])[k], population.row(parent1)[296])), 1/population.row(parent1)[296]);
				}
				if(betterMatch < bestMatch)
				{
					bestMatch = betterMatch;
					bestIndex = j;
				}
			}
			
			//Mate
			for(int j = 0; j < population.cols(); j++)
			{
				double pickOne = r.nextDouble();
				if(pickOne < 0.5)
					population.row(dead[i])[j] = population.row(parent1)[j];
				else
					population.row(dead[i])[j] = population.row(parent2[bestIndex])[j];
			}
		}
	}
	
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

	
	// ----------------------------------------------------------------
	// The contents of this file are distributed under the CC0 license.
	// See http://creativecommons.org/publicdomain/zero/1.0/
	// ----------------------------------------------------------------


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


	static void fullTournament(double [] w) throws Exception {
		ArrayList<IAgent> al = new ArrayList<IAgent>();
		//al.add(new PrescientMoron());
		//al.add(new Mixed());
		//al.add(new Human());
		//al.add(new Blitz());
		//al.add(new SittingDuck());
		//al.add(new AggressivePack());
		al.add(new Winner2015a());
		al.add(new Winner2015b());
		al.add(new Winner2016a());
		al.add(new BellBrent());
		al.add(new BellBrent2(w));

		Controller.doTournament(al);
	}
	
	static double[] evolveWeights() throws Exception
	{
		// Create a random initial population
		Random r = new Random();
		Matrix population = new Matrix(100, 297);
		for(int i = 0; i < 100; i++)
		{
			double[] chromosome = population.row(i);
			for(int j = 0; j < chromosome.length; j++)
				chromosome[j] = r.nextGaussian();
		}
		System.out.println(Arrays.toString(population.row(0)));
		// Evolve the population
		
		//META Parameters
		double p = 0.9; //Probability of mutation
		double d = 0.08; //Std Deviation for mutation
		double survival = 0.8; //Probability winner survives (Natural Selection)
		int numberDead = 5; //Number that die from battle
		int numberOfSecondParents = 20; //Number of second parents for replinishment
		int powerDifference = 2; //Power difference when finding distance between parents
		
		for(int i = 0; i < 100; i++)
		{
			population.row(i)[291] = p; //Probability of mutation
			population.row(i)[292] = d;
			population.row(i)[293] = survival;
			population.row(i)[294] = numberDead;
			population.row(i)[295] = numberOfSecondParents;
			population.row(i)[296] = powerDifference;
		}
		
		//Evolve for a long time
		for(int repeat = 0; repeat < 100; repeat++)
		{
			//Reduce probability of mutation slowly
			for(int i = 0; i < population.rows(); i++)
			{
				if(population.row(i)[291] > 0.1 && population.row(i)[292] > 0.1)
				{
					population.row(i)[291] -= 0.6/100;
					population.row(i)[292] -= 0.9/100;
				}
			}
			if((double)repeat / 1 % 5 == 0)
			{
				System.out.println("Percent Complete: " + (double)repeat / 1);
			}
			
			//Mutate
			//System.out.println("Before Mutate: " + Arrays.toString(population.row(0)));
			Mutate(population);
			//System.out.println("After Mutate:  " + Arrays.toString(population.row(0)));
			
			//Natural Selection
			int[] dead = NaturalSelection(population);
			//for(int i = 0; i < dead.length; i++)
			//{
			//	if (dead[i] == 0)
			//			System.out.println("Population 0 Died");
			//}

			//Replinish
			//System.out.println("Dead: " + Arrays.toString(dead));
			//System.out.println("Before Replenish: " + Arrays.toString(population.row(dead[0])));
			Replenish(population, dead);
			//System.out.println("After Replenish:  " + Arrays.toString(population.row(dead[0])));
			//System.out.println("After Replenish:  " + Arrays.toString(population.row(0)));

		}
		//System.out.println("");
		//System.out.println(population.toString());

		// Return an arbitrary member from the population
		double[][] returnDoubles = new double[100][291];
		for(int j = 0; j < 100; j++)
		{
			for(int i = 0; i < 291; i++)
			{
				returnDoubles[j][i] = population.row(j)[i];
			}
		}
		int whichWon = 0;
		int thisoneWon = 0;
		for(int i = 0; i < 100; i++)
		{
			try {
				int chooseOpponent = r.nextInt(4);
				IAgent agent;
				if(chooseOpponent == 0)
					agent = new Winner2016a();
				else if(chooseOpponent == 1)
					agent = new Winner2015a();
				else if(chooseOpponent == 1)
					agent = new Winner2015b();
				else
					agent = new BellBrent();
				if(r.nextDouble() < 0.5)
					whichWon = Controller.doBattleNoGui(new BellBrent2(), agent);
				else
					whichWon = Controller.doBattleNoGui(agent, new BellBrent2());
			} catch (Exception e)
			{
				System.out.println("Exception: " + e.getMessage());
			}
			if (whichWon == -1)
			{
				//System.out.println("Red Won!!");
				thisoneWon = i;
			}
				
		}
		return returnDoubles[thisoneWon];
	}

	public static void main(String[] args) throws Exception {
		//Controller.doBattle(new BellBrent(), new Winner2016a());
		//Controller.doBattle(new Mixed(), new Blitz());
		//Controller.doBattle(new Mixed(), new AggressivePack());
		//Controller.doBattle(new Blitz(), new Mixed());
		//Controller.doBattle(new Human(), new SittingDuck());
		//Controller.doBattle(new Mixed(), new SittingDuck());
		//Controller.doBattle(new Human(), new Blitz());
		//Controller.doBattle(new BellBrent(), new SittingDuck());
		//Controller.doBattle(new PrescientMoron(), new Human());
		//Controller.doBattle(new BellBrent(), new PrescientMoron());
		//fullTournament();
		//double[] w = evolveWeights();
		//System.out.println(w);
		Controller.doBattle(new BellBrent2(), new Winner2016a());
		//fullTournament(w);

	}
}
