package thud.simpel;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.TreeMap;

public class SimpelModel
{	
	private String directory = "/home/kerryn/git/2023-06-SimpelSoilWaterModels/SIMPLE/";
	
	public final static int PRECIPITATION=0; 
	public final static int SNOW_WATER_EQUI=1; 
	public final static int SNOW_MELT_RAIN=2;
    public final static int ETP_COEFF=3; 
    public final static int ETP_INPUT=4; 
    public final static int LAI=5; 
    public final static int I_CAP=6; 
    public final static int INT_ETI_LEAF=7;
    public final static int I_BAL=8;
    public final static int I_PREC=9;
    public final static int I_REM=10; 
    public final static int I_ETI_LITTER=11; 
    public final static int BILANZ=12; 
    public final static int CONTENT=13; 
    public final static int S_RESTN=14; 
    public final static int INF_LIMIT=15; 
    public final static int P_LINF=16; 
    public final static int REST_ETA=17; 
    public final static int BALANCE_SOIL=18;
    public final static int ETA=19; 
    public final static int ET_BALANCE=20; 
    public final static int SEEPAGE=21;	
    public final static int STORAGE=22;	
    public final static int SURFACE_RUNOFF=23;
    public final static int RUNOFF_TOTAL=24; 
    public final static int I_LEAF=25;	
    public final static int I_LITTER=26;	
    public final static int ETA_TOTAL=27;	
    
    public final static int INPUT_DATE=0;
    public final static int INPUT_DOY=1;	
    public final static int INPUT_P=2;	
    public final static int INPUT_T14=3;	
    public final static int INPUT_R14=4;	
    public final static int INPUT_ET0=5;
    
    public final static int LAI_DOY = 0;
    public final static int LAI_LAI = 1;
    public final static int LAI_V3 = 2;
    public final static double SECONDS_PER_DAY = 86400.;
    
    public final String tab = "\t";
    public final String linefeed = "\n";
    private ArrayList<String> bucket_model_dates;
    
  //public int numColumnsInput = 6;  //hardcoding this for now.
    private boolean readETfromFile = false; //hardcoding this for now. $kf2023-06-23
    private double ts = 86400.; // model time step in seconds
    int krowEnd = 0;

	public static void main(String[] args)
	{
		SimpelModel s = new SimpelModel();
		if (args.length>0)
		{
			String arg0 = args[0];
			s.setDirectory(arg0);
		}
		if (args.length>1)
		{
			String arg1 = args[1];
			s.setReadETFromFile(arg1);
		}
		s.process();
	}
	
	public void process()
	{
		ArrayList<String[]> Input = readLargerTextFileAlternateToArray(directory + "Global_Input.txt", true, "\\s+");
		ArrayList<String[]> Landuse = readLargerTextFileAlternateToArray(directory + "Land_Use.txt", true, "\\s+");
		ArrayList<String[]> Landuses = readLargerTextFileAlternateToArray(directory + "Land_Use.txt", false, "\\s+");
		String[] landusesHeader = Landuses.get(0);
		TreeMap<String,Integer> landusesIndexes = new TreeMap<String,Integer>();
		for (int i=1;i<landusesHeader.length;i++)
		{
			landusesIndexes.put( landusesHeader[i], i);
		}
		ArrayList<String[]> LAI_model = readLargerTextFileAlternateToArray(directory + "LAI_model.txt", true, "\\s+");
		ArrayList<String[]> Soil = readLargerTextFileAlternateToArray(directory + "Soil_physics.txt", false, ";");
		String output = directory + "SIMPLE_bucket_model.csv";

		double[][] bucket_model = SIMPLE_function(Input, Landuse, LAI_model, Soil, landusesIndexes);
		write_csv(bucket_model, output);

	}
	
	public double[][] SIMPLE_function(ArrayList<String[]> Input, ArrayList<String[]> Landuse, 
			ArrayList<String[]> LAI_model, ArrayList<String[]> Soil, TreeMap<String,Integer> landusesIndexes)
	{		
		double Field_Capacity_Percent=Double.parseDouble(Soil.get(0)[1]);//Soil[1][2]
		double Permanent_Wilting_Point=Double.parseDouble(Soil.get(1)[1]);//Soil[2][2]
		double Start_of_Reduction_Percent=Double.parseDouble(Soil.get(2)[1]);//Soil[3][2]
		double Root_Depth=Double.parseDouble(Soil.get(3)[1]);//Soil[4][2]
		double Init_Value_Soil_Percent=Double.parseDouble(Soil.get(4)[1]);//Soil[5][2]
		double Field_Capacity=Double.parseDouble(Soil.get(5)[1]);//Soil[6][2]
		double Perm_Wilting_Point=Double.parseDouble(Soil.get(6)[1]);//Soil[7][2]
		double FWC=Double.parseDouble(Soil.get(7)[1]);//Soil[8][2]
		double Start_of_Reduction=Double.parseDouble(Soil.get(8)[1]);//Soil[9][2]
		double Init_Value_Soil=Double.parseDouble(Soil.get(9)[1]);//Soil[10][2]
		double Depth_of_soil=Double.parseDouble(Soil.get(10)[1]);//Soil[11][2]
		String Land_use=Soil.get(11)[1];//Soil[12][2]
		double Minimum_LAI=Double.parseDouble(Soil.get(12)[1]);//Soil[13][2]
		double Maximum_LAI=Double.parseDouble(Soil.get(13)[1]);//Soil[14][2]
		double Direct_runoff_factor=Double.parseDouble(Soil.get(14)[1]);//Soil[15][2]			
		double Koeff_c=Double.parseDouble(Soil.get(15)[1]);//Soil[16][2]
		double Lambda=Double.parseDouble(Soil.get(16)[1]);//Soil[17][2]
		double Cap_Litter=Double.parseDouble(Soil.get(17)[1]);//Soil[18][2]
		double Init_Value_Litter=Double.parseDouble(Soil.get(18)[1]);//Soil[19][2]
		double Litter_Reduction_factor=Double.parseDouble(Soil.get(19)[1]);//Soil[20][2]
		
		int landuseIndex = landusesIndexes.get(Land_use); 
		
//		# SIMPLE_function
//		# original Excel spreadsheet developed by Georg Hörmann
//		# extended (snow, surface runoff) by Kristian Förster 2022
//		# translated to R by Zoe Bovermann 2023
//      # translated to Java by Kerry Nice 20 June 2023

//		  # prepare bucket model matrix
		double[][] bucket_model = new double[Input.size()][28];
		bucket_model_dates = new ArrayList<String>();
		
//		  # copy Input into bucket-model
//		  #col 1: A Date
//		  # col 2: B Precipitation				  
		  for (int i=0;i<Input.size();i++)
		  {
			  bucket_model_dates.add(Input.get(i)[INPUT_DATE]);
			  bucket_model[i][PRECIPITATION]=Double.parseDouble(Input.get(i)[INPUT_P]);
		  }
		  
//		  # water balance check
		  double sum_prec   = 0.;
		  double sum_etr    = 0.;
		  double sum_runoff = 0.;
		  double init_stor  = Init_Value_Soil; 
		  double init_swe   = 0.;
		    
		  double[][] laiModel = new double[4][3];
		  for (int i=0;i<LAI_model.size();i++)
		  {
			  laiModel[i][LAI_DOY] = Double.parseDouble(LAI_model.get(i)[0]);
			  laiModel[i][LAI_LAI] = Double.parseDouble(LAI_model.get(i)[1]);
		  }

//		  # update LAI model according to soil physics ($kf 2023-03-28) -> ignore land use table
//		  # update Litter according to land use table
		  laiModel[0][LAI_V3]= Double.parseDouble(Soil.get(12)[1]);
		  laiModel[1][LAI_V3]= Double.parseDouble(Soil.get(13)[1]);
		  laiModel[2][LAI_V3]= Double.parseDouble(Soil.get(13)[1]);
		  laiModel[3][LAI_V3]= Double.parseDouble(Soil.get(12)[1]);
//		  # update Litter according to land use table
		  double landuse_DayDegree = Double.parseDouble(Landuse.get(16)[landuseIndex]);
		  
		  double laimodel13 = Double.parseDouble(LAI_model.get(0)[2]); // LAI_model[1,3]
		  double laimodel12 = Double.parseDouble(LAI_model.get(0)[1]); // LAI_model[1,2]
		  
		  double laimodel23 = Double.parseDouble(LAI_model.get(1)[2]); // LAI_model[2,3]
		  double laimodel22 = Double.parseDouble(LAI_model.get(1)[1]); // LAI_model[2,2]
		  
		  double laimodel33 = Double.parseDouble(LAI_model.get(2)[2]); // LAI_model[3,3]
		  double laimodel32 = Double.parseDouble(LAI_model.get(2)[1]); // LAI_model[3,2]
		  
		  double laimodel43 = Double.parseDouble(LAI_model.get(3)[2]); // LAI_model[4,3]
		  double laimodel42 = Double.parseDouble(LAI_model.get(3)[1]); // LAI_model[4,2]
        
          double nt = ts / SECONDS_PER_DAY; // fractional time step
		   
//		  # calculate the rest
		  for(int krow=0;krow<Input.size();krow++)
		  {			  
			  double t14 = Double.parseDouble(Input.get(krow)[INPUT_T14]);
			  double rh14 = Double.parseDouble(Input.get(krow)[INPUT_R14]); // $kf: added for internal ETP calculation
			  double doy = Double.parseDouble(Input.get(krow)[INPUT_DOY]);// Input[krow][INPUT_DOY]
			  double p = Double.parseDouble(Input.get(krow)[INPUT_P]);// Input[krow][INPUT_P]
//		    # col 3+4:  
//		    # first row
		    if(krow==0)
		    {    	
//		      # col 3: C snow water equi
		    	bucket_model[krow][SNOW_WATER_EQUI] = init_swe; //# col 3 // $kf set to init_swe for consistency with water balance
		      
//		      # col 4: D Snow melt + rain
		      double[] temp = new double[]{0, landuse_DayDegree*t14 * nt};
		      if(t14 >= 0)
		      {
		        bucket_model[krow][SNOW_MELT_RAIN] = Math.min(0,(max(temp))) + bucket_model[krow][PRECIPITATION];
		      }
		      else
		      {
		        bucket_model[krow][SNOW_MELT_RAIN] = Math.min(0,(max(temp))) + 0;
		      }
		      
//		      # col 3 (dependent on col 4)
		      bucket_model[krow][SNOW_WATER_EQUI] = bucket_model[krow][SNOW_WATER_EQUI] + bucket_model[krow][PRECIPITATION] 
		    		  - bucket_model[krow][SNOW_MELT_RAIN];
		    }
		    
//		    # rest of the rows
		    if(krow > 0) // $kf: applies to all rows greater than zero
		    {
		    	double t14_1 = Double.parseDouble(Input.get(krow)[INPUT_T14]);
		    	
		      double[] temp = new double[]{0, landuse_DayDegree*t14_1 * nt};
		      if(t14_1 < 0)
		      {
//		        # col 4: Snow melt + rain 
		        bucket_model[krow][SNOW_MELT_RAIN] = 
		        		Math.min(max(temp), bucket_model[krow-1][SNOW_WATER_EQUI]+bucket_model[krow][PRECIPITATION]);
		        
//		        # col 3 (dependent on col 4)
		        bucket_model[krow][SNOW_WATER_EQUI] = bucket_model[krow-1][SNOW_WATER_EQUI] + bucket_model[krow][PRECIPITATION] 
		        		- bucket_model[krow][SNOW_MELT_RAIN]+0;
		      }
		      else
		      {
		        bucket_model[krow][SNOW_MELT_RAIN] = Math.min(max(temp),(bucket_model[krow-1][SNOW_WATER_EQUI])+0) 
		        		+ bucket_model[krow][PRECIPITATION];
		        bucket_model[krow][SNOW_WATER_EQUI] = bucket_model[krow-1][SNOW_WATER_EQUI] + 0 
		        		- bucket_model[krow][SNOW_MELT_RAIN] +  bucket_model[krow][PRECIPITATION];
		      }
		    }
//		    # col 5: E ETP Coeff.Landuse[17,which(colnames(Landuse)==Soil[12,2])]
		    String date = bucket_model_dates.get(krow);
		    String[] dateSplit = date.split("\\.");
		    int month = Integer.parseInt(dateSplit[1]);
		    bucket_model[krow][ETP_COEFF] = Double.parseDouble(Landuse.get(month-1)[landuseIndex]);
		    
//		    # col 6: F ETP Input
		    if(readETfromFile) // if(numColumnsInput >= 6) // $kf
		    {   
		    	  //# check if input variable is given
//		          # if yes, use given ETP0
		    	 double et0_1 = Double.parseDouble(Input.get(krow)[INPUT_ET0]);
			      bucket_model[krow][ETP_INPUT] = et0_1;
			      if(krow==2)
			      {
			    	  System.out.println("Reading external ETP from file...");
			      }
		    }
		    else 
		    { 
		    	//# if not calculate ETP
		    	bucket_model[krow][ETP_INPUT] = Math.min(7, bucket_model[krow][ETP_COEFF]*6.11*
			    		  Math.pow(10,((7.5*t14)/(237.3+t14)))*(1.-(rh14/100.))) * nt; // $kf: relative humidity fixed
		    }
		    
//		    # col 7: G LAI
		    if(doy <= laimodel12)
		    {
		      bucket_model[krow][LAI] = laimodel13;
		    }
		    else if(doy <= laimodel22)
		    {
		      bucket_model[krow][LAI] = Minimum_LAI+(laimodel23-laimodel13)*((doy-laimodel12)/(laimodel22-laimodel12));
		    }
		    else if(doy <= laimodel32)
		    {
		      bucket_model[krow][LAI] = laimodel33;
		    }
		    else if(doy <= laimodel42)
		    {
		      bucket_model[krow][LAI] = Maximum_LAI+(laimodel43-laimodel33)*((doy-laimodel32)/(laimodel42-laimodel32));
		    }
		    else
		    {
		      bucket_model[krow][LAI] = laimodel43;
		    }

//		    # col 8: H I-Cap
		    bucket_model[krow][I_CAP] = 0.35*bucket_model[krow][LAI];

//		    # col 9: I ETi 
		    bucket_model[krow][INT_ETI_LEAF] = Math.min(bucket_model[krow][I_CAP], bucket_model[krow][ETP_INPUT]);

//		    # col 10: J I-Bal.   
		    bucket_model[krow][I_BAL] = bucket_model[krow][SNOW_MELT_RAIN] - bucket_model[krow][INT_ETI_LEAF];
		    
//		    # col 11: K I-Prec.   
		    bucket_model[krow][I_PREC] = Math.max(bucket_model[krow][I_BAL],0) ;
		    
//		    # col 12: L I-Rem.
		    bucket_model[krow][I_REM] = -1*Math.min(bucket_model[krow][I_BAL],0)
		    		-bucket_model[krow][INT_ETI_LEAF]+bucket_model[krow][ETP_INPUT];

//		    # col 13+14+15
//		    # first row
		    if(krow==0)
		    {
//		      # col 13: M ETi Litter
		      double[] temp = new double[]{Cap_Litter,(0+bucket_model[krow][I_PREC])/Litter_Reduction_factor};
		      bucket_model[krow][I_ETI_LITTER] = Math.min(bucket_model[krow][I_REM],min(temp));
		      
//		      # col 14: N Bilanz (needed for col 13 & 15)
		      bucket_model[krow][BILANZ] = 0 + bucket_model[krow][I_PREC] - bucket_model[krow][I_ETI_LITTER];
		      
//		      # col 15: O Content (needed for col 13 & 14)
		      if(bucket_model[krow][BILANZ] > Cap_Litter)
		      {
		        bucket_model[krow][CONTENT] = Cap_Litter;
		      }
		      else 
		      {
		    	  bucket_model[krow][CONTENT] = Math.max(0,bucket_model[krow][BILANZ]);
		      }		      
		    }
		    else
		    {//# rest of the rows
//		      # col 13: M ETi Litter
		      double[] temp = new double[]{Cap_Litter, (bucket_model[krow-1][CONTENT]+bucket_model[krow][I_PREC])/Litter_Reduction_factor};
		      bucket_model[krow][I_ETI_LITTER] = Math.min(bucket_model[krow][I_REM],min(temp));
		      
//		      # col 14: N Bilanz (needed for col 13 & 15)
		      bucket_model[krow][BILANZ] = bucket_model[krow-1][CONTENT]  + bucket_model[krow][I_PREC] - bucket_model[krow][I_ETI_LITTER];
		      
//		      # col 15: O Content (needed for col 13 & 14)
		      if(bucket_model[krow][BILANZ] > Cap_Litter)
		      {
		        bucket_model[krow][CONTENT] = Cap_Litter;
		      }
		      else
		      {
		    	  bucket_model[krow][CONTENT] = Math.max(0,bucket_model[krow][BILANZ]);
		      }
		    }		    
//		    # col 16: P S-REstn
		    bucket_model[krow][S_RESTN] = Math.max(0, (bucket_model[krow][BILANZ]-bucket_model[krow][CONTENT]));
		    
//		    # col 17-24:
//		    # first row
		    if(krow==0)
		    {
//		      # col 17: Q Inf-Limit
		      bucket_model[krow][INF_LIMIT] = (Field_Capacity-Init_Value_Soil)*0.25*(1.-Direct_runoff_factor/100.);
		      
//		      # col 18: R P-Inf
		      bucket_model[krow][P_LINF] = Math.min(bucket_model[krow][S_RESTN],bucket_model[krow][INF_LIMIT])*(1.-Direct_runoff_factor/100.);
		      
//		      # col 19: S S-Rest
		      bucket_model[krow][REST_ETA] = -1*Math.min(0,bucket_model[krow][BILANZ])+bucket_model[krow][I_REM]-bucket_model[krow][I_ETI_LITTER];
		      
//		      # col 20: T Balance soil
		      bucket_model[krow][BALANCE_SOIL] = Init_Value_Soil + bucket_model[krow][P_LINF];
		      
//		      # col 21: U ETa
		      if(bucket_model[krow][BALANCE_SOIL] > Start_of_Reduction)
		      {
		        bucket_model[krow][ETA] =  bucket_model[krow][REST_ETA];
		      }
		      else
		      {
		        bucket_model[krow][ETA] =  bucket_model[krow][REST_ETA]*(bucket_model[krow][BALANCE_SOIL]-Perm_Wilting_Point)/
		          (Start_of_Reduction- Perm_Wilting_Point);
		      }		      
//		      # col 22: V ET-Balance
		      bucket_model[krow][ET_BALANCE] = bucket_model[krow][BALANCE_SOIL] - bucket_model[krow][ETA];
		      
//		      # col 23: W Seepage
		      if(bucket_model[krow][ET_BALANCE] <= Field_Capacity)
		      {
		        bucket_model[krow][SEEPAGE] =  Lambda*Math.pow( (bucket_model[krow][ET_BALANCE]-Perm_Wilting_Point),2);
		      }
		      else
		      {
		        bucket_model[krow][SEEPAGE] =  Lambda*Math.pow( (Field_Capacity-Perm_Wilting_Point),2);
		      }		      
//		      # col 24: X Storage Init.-Value
		      if(bucket_model[krow][ET_BALANCE] > Field_Capacity)
		      {
		        bucket_model[krow][STORAGE] =  Field_Capacity;
		      }
		      else
		      {
		        bucket_model[krow][STORAGE] = bucket_model[krow][ET_BALANCE] - bucket_model[krow][SEEPAGE];
		      }		      
		    }
		    else
		    {
//		      # col 17: Q Inf-Limit
		      bucket_model[krow][INF_LIMIT] = (Field_Capacity-bucket_model[krow-1][STORAGE])*0.25;
		      
//		      # col 18: R P-Inf
		      bucket_model[krow][P_LINF] = Math.min(bucket_model[krow][S_RESTN],bucket_model[krow][INF_LIMIT])*(1.-Direct_runoff_factor/100.);
		      
//		      # col 19: S S-Rest
		      bucket_model[krow][REST_ETA] = -1*Math.min(0,bucket_model[krow][BILANZ])+bucket_model[krow][I_REM]-bucket_model[krow][I_ETI_LITTER];
		      
//		      # col 20: T Balance soil
		      bucket_model[krow][BALANCE_SOIL] = bucket_model[krow-1][STORAGE] + bucket_model[krow][P_LINF];
		      
//		      # col 21: U ETa
		      if(bucket_model[krow][BALANCE_SOIL] > Start_of_Reduction)
		      {
		        bucket_model[krow][ETA] =  bucket_model[krow][REST_ETA];
		      }
		      else
		      {
		        bucket_model[krow][ETA] =  bucket_model[krow][REST_ETA]*(bucket_model[krow][BALANCE_SOIL]-Perm_Wilting_Point)/
		          (Start_of_Reduction- Perm_Wilting_Point);
		      }
		      
//		      # col 22: V ET-Balance
		      bucket_model[krow][ET_BALANCE] = bucket_model[krow][BALANCE_SOIL] - bucket_model[krow][ETA];
		      
//		      # col 23: W Seepage
		      if(bucket_model[krow][ET_BALANCE] <= Field_Capacity)
		      {
		        bucket_model[krow][SEEPAGE] =  Lambda*
		        		Math.pow((bucket_model[krow][ET_BALANCE]-Perm_Wilting_Point),2);
		      }
		      else
		      {
		        bucket_model[krow][SEEPAGE] =  Lambda*Math.pow((Field_Capacity-Perm_Wilting_Point),2);
		      }
		      
//		      # col 24: X Storage Init.-Value
		      if(bucket_model[krow][ET_BALANCE] > Field_Capacity)
		      {
		        bucket_model[krow][STORAGE] =  Field_Capacity;
		      }
		      else
		      {
		        bucket_model[krow][STORAGE] = bucket_model[krow][ET_BALANCE] - bucket_model[krow][SEEPAGE];
		      }
		    }
//		    # col 25: Y surface runoff
		    if(bucket_model[krow][ET_BALANCE] > Field_Capacity)
		    {
		      bucket_model[krow][SURFACE_RUNOFF] = bucket_model[krow][ET_BALANCE]-Field_Capacity+bucket_model[krow][S_RESTN]
		    		  -bucket_model[krow][P_LINF];
		    }
		    else
		    {
		        bucket_model[krow][SURFACE_RUNOFF] = bucket_model[krow][S_RESTN]-bucket_model[krow][P_LINF];
		    }
		    
//		    # col 26: Z Runofftotal
		    bucket_model[krow][RUNOFF_TOTAL] = bucket_model[krow][SURFACE_RUNOFF] + bucket_model[krow][SEEPAGE];
		    
//		    # col 27: AA I-Leaf
		    bucket_model[krow][I_LEAF] = bucket_model[krow][ETP_INPUT] -  bucket_model[krow][I_REM];
		    
//		    # col 28: AB I-Litter
		    bucket_model[krow][I_LITTER] = bucket_model[krow][I_REM] -  bucket_model[krow][REST_ETA];
		    
//		    # col 29: AC ETa Total
		    bucket_model[krow][ETA_TOTAL] = bucket_model[krow][I_LITTER] + bucket_model[krow][I_LEAF] + bucket_model[krow][ETA];

		    sum_prec = sum_prec + p;
		    sum_etr  = sum_etr  + bucket_model[krow][ETA_TOTAL];
		    sum_runoff = sum_runoff + bucket_model[krow][RUNOFF_TOTAL];
		    krowEnd = krow;
		  }
//		  # water balance check
		  double water_balance = sum_prec - sum_etr - sum_runoff + init_swe + init_stor - bucket_model[krowEnd][STORAGE] - bucket_model[krowEnd][SNOW_WATER_EQUI];
		  System.out.println("water balance check");
		  System.out.println(water_balance);
		  
		  return bucket_model;
		}
	
	double max(double[] items)
	{
		return Math.max(items[0], items[1]);
	}
	double min(double[] items)
	{
		return Math.min(items[0], items[1]);
	}
	
   final static Charset ENCODING = StandardCharsets.UTF_8;	
   public ArrayList<String[]> readLargerTextFileAlternateToArray(String aFileName, boolean skipHeader, String deliminator) 
   {
	   ArrayList<String[]> returnFile = new ArrayList<String[]>();
	   int count = 0;
	    Path path = Paths.get(aFileName);
	    try (BufferedReader reader = Files.newBufferedReader(path, ENCODING))
	    {
	      String line = null;
	      while ((line = reader.readLine()) != null) 
	      {
	    	  if (skipHeader && count==0)
	    	  {
	    		  count++;
	    		  continue;
	    	  }
	    	  if (line != null && line.trim().equals(""))
	    	  {
	    		  continue;
	    	  }
	    	  String[] lineSplit = line.split(deliminator);
	    	  returnFile.add(lineSplit);
	    	  count ++;
	      }      
	    }
		catch (IOException e)
		{			
			e.printStackTrace();
		}
	    return returnFile;
	  }
	public void writeFile(String text, String filename)
	{	
		FileOutputStream out; 
		PrintStream p; 
		try
		{
			out = new FileOutputStream(filename);
			p = new PrintStream(out);
			p.println(text);
			p.close();
		} catch (Exception e)
		{
			System.err.println("Error writing to file");
		}
		p=null;
		out=null;
	}	   
	public void write_csv(double[][] bucket_model, String output)
	{
//		Date	Precipitation	Snow_Water_Equi	Snow_melt+rain	ETP_Coeff	ETP_Input	LAI	I-Cap	Int. ETi(leaf)	I-Bal	I-Prec	I-Rem	I_ETi(litter)	Bilanz	Content	S-REstn	Inf-Limit	P-linf	Rest-ETA	Balance_soil	ETa	ET-Balance	Seepage	Storage	Surface_runoff	Runoff_total	I-Leaf	I-Litter	ETaTotal
		String header = "Date" + tab
				+ "Precipitation" + tab
				+ "Snow_Water_Equi" + tab
				+ "Snow_melt+rain" + tab
				+ "ETP_Coeff" + tab
				+ "ETP_Input" + tab
				+ "LAI" + tab
				+ "I-Cap" + tab
				+ "Int. ETi(leaf)" + tab
				+ "I-Bal" + tab
				+ "I-Prec" + tab
				+ "I-Rem" + tab
				+ "I_ETi(litter)" + tab
				+ "Bilanz" + tab
				+ "Content" + tab
				+ "S-REstn" + tab
				+ "Inf-Limit" + tab
				+ "P-linf" + tab
				+ "Rest-ETA" + tab
				+ "Balance_soil" + tab
				+ "ETa" + tab
				+ "ET-Balance" + tab
				+ "Seepage" + tab
				+ "Storage" + tab
				+ "Surface_runoff" + tab
				+ "Runoff_total" + tab
				+ "I-Leaf" + tab
				+ "I-Litter" + tab
				+ "ETaTotal";
		StringBuilder sb = new StringBuilder();
		sb.append(header + linefeed);
		int width = bucket_model[0].length;
		int height = bucket_model.length;
		for (int i=0;i<height;i++)
		{
			sb.append(bucket_model_dates.get(i));
			for (int j=0;j<width;j++)
			{
				sb.append(tab + bucket_model[i][j]);
			}
			sb.append(linefeed);
		}
		writeFile(sb.toString(), output);
	}

	public String getDirectory()
	{
		return directory;
	}
	public void setDirectory(String directory)
	{
		this.directory = directory;
	}	
	public boolean getReadETFromFile()
	{
		return readETfromFile;
	}
	public void setReadETFromFile(String readET)
	{
		if (readET.equals("readET"))
		{
			readETfromFile=true;
		}
		else
		{
			readETfromFile=false;
		}
	}
}
