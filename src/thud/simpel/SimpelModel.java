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
	  //TODO utilize kdown, ldown, wind speed, now available in Global_Input
	  //TODO, Use timestep values, Hourly value for timestep now configured in Soil_physics: Timestep;24;
	
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
    
    public final static int K_DOWN=28;	
    public final static int WIND_SPEED=29;	
    public final static int L_DOWN=30;	
    
    public static int INPUT_DATE=0;
    public static int INPUT_DOY=1;	
    public static int INPUT_P=2;	
    public static int INPUT_T14=3;	
    public static int INPUT_R14=4;	
    public static int INPUT_ET0=5;    
    public static int INPUT_K_DOWN=6;
    public static int INPUT_WIND_SPEED=7;
    public static int INPUT_L_DOWN=8;
    
    public final static int LAI_IND = 0;
    public final static int LAI_DOY = 1;
    public final static int LAI_PHASE = 2;
    public final static int LAI_LAI = 3;
    public final static double SECONDS_PER_DAY = 86400.;
    
    
    
    public final static String FIELD_CAPACITY_PERCENT = "Field Capacity %";
    public final static String PERMANENT_WILTING_POINT = "Permanent Wilting Point %";
    public final static String START_OF_REDUCTION_PERCENT = "Start of Reduction %";
    public final static String ROOT_DEPTH = "Root Depth";
    public final static String INIT_VALUE_SOIL_PERCENT = "Init-Value Soil %";
    //public final static String FIELD_CAPACITY = "Field Capacity";
    //public final static String PERM_WILTING_POINT = "Perm. Wilting Point";
    //public final static String FWC = "FWC";
    //public final static String START_OF_REDUCTION = "Start of Reduction";
    //public final static String INIT_VALUE_SOIL = "Init-Value Soil";
    //public final static String DEPTH_OF_SOIL = "Depth of soil";	
    public final static String LAND_USE = "Land use";
    public final static String MINIMUM_LAI = "Minimum LAI";
    public final static String MAXIMUM_LAI = "Maximum LAI";
    public final static String INTC_COVERED_FRACTION = "Vegetation Fraction";
    public final static String INTC_LAYER_THICKNESS = "Layer Thickness";
    public final static String INTC_DRAINAGE_EXP_B = "Drainage Coeff. b";
    public final static String INTC_DRAINAGE_MAX = "Max. Drainage Rate";
    public final static String DIRECT_RUNOFF_FACTOR = "Direct runoff factor";				
    public final static String GLUGLA_C = "Glugla coeff.";
    //public final static String LAMBDA = "Lambda";	
    public final static String CAP_LITTER = "Cap. Litter";
    public final static String INIT_VALUE_LITTER = "Init-Value Litter";
    public final static String LITTER_REDUCTION_FACTOR = "Litter Reduction factor";	
    public final static String TIMESTEP = "Timestep";	
    
    public final String tab = "\t";
    public final String linefeed = "\n";
    private ArrayList<String> bucket_model_dates;
    
  //public int numColumnsInput = 6;  //hardcoding this for now.
    private boolean readETfromFile = false; //hardcoding this for now. $kf2023-06-23
    // private double ts = 86400.; // model time step in seconds
    int krowEnd = 0;
    boolean hasKdown = false;
    boolean hasWindSpeed = false;
    boolean hasLdown = false;

    //public double INTC_DRAINAGE_MAX_D = 2.88; // mm/d
    //public double INTC_DRAINAGE_EXP_B = 3.7;  // mm^(-1)

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
		ArrayList<String[]> Input = readLargerTextFileAlternateToArray(directory + "Global_Input.txt", true, "\\s+",false);
		ArrayList<String[]> InputHeader = readLargerTextFileAlternateToArray(directory + "Global_Input.txt", false, "\\s+",true);
		matchInputHeaders(InputHeader);
		ArrayList<String[]> Landuse = readLargerTextFileAlternateToArray(directory + "Land_Use.txt", true, "\\s+",false);
		ArrayList<String[]> Landuses = readLargerTextFileAlternateToArray(directory + "Land_Use.txt", false, "\\s+",false);
		String[] landusesHeader = Landuses.get(0);
		TreeMap<String,Integer> landusesIndexes = new TreeMap<String,Integer>();
		for (int i=1;i<landusesHeader.length;i++)
		{
			landusesIndexes.put( landusesHeader[i], i);
		}
		ArrayList<String[]> LAI_model = readLargerTextFileAlternateToArray(directory + "LAI_model.txt", true, "\\s+",false);
		TreeMap<String,String[]> Soil = readLargerTextFileAlternateToTreeMap(directory + "Soil_physics.txt", false, ";");
		String output = directory + "SIMPLE_bucket_model.csv";

		double[][] bucket_model = SIMPLE_function(Input, Landuse, LAI_model, Soil, landusesIndexes);
		write_csv(bucket_model, output);

	}
	
	
	public double[][] SIMPLE_function(ArrayList<String[]> Input, ArrayList<String[]> Landuse, 
			ArrayList<String[]> LAI_model, TreeMap<String,String[]> Soil, TreeMap<String,Integer> landusesIndexes)
	{		
		double Field_Capacity_Percent=Double.parseDouble(Soil.get(FIELD_CAPACITY_PERCENT)[1]);//Soil[1][2]
		double Permanent_Wilting_Point=Double.parseDouble(Soil.get(PERMANENT_WILTING_POINT)[1]);//Soil[2][2]
		double Start_of_Reduction_Percent=Double.parseDouble(Soil.get(START_OF_REDUCTION_PERCENT)[1]);//Soil[3][2]
		double Root_Depth=Double.parseDouble(Soil.get(ROOT_DEPTH)[1]);//Soil[4][2]
		double Init_Value_Soil_Percent=Double.parseDouble(Soil.get(INIT_VALUE_SOIL_PERCENT)[1]);//Soil[5][2]
		double Field_Capacity= Field_Capacity_Percent * Root_Depth / 10.;//Double.parseDouble(Soil.get(FIELD_CAPACITY)[1]);//Soil[6][2]
		double Perm_Wilting_Point=Permanent_Wilting_Point * Root_Depth / 10.;//Double.parseDouble(Soil.get(PERM_WILTING_POINT)[1]);//Soil[7][2]
		double fwc=Field_Capacity-Perm_Wilting_Point;//Double.parseDouble(Soil.get(FWC)[1]);//Soil[8][2]
		double Start_of_Reduction= Start_of_Reduction_Percent * Root_Depth / 10.;///Double.parseDouble(Soil.get(START_OF_REDUCTION)[1]);//Soil[9][2]
		double Init_Value_Soil=Init_Value_Soil_Percent * Root_Depth / 10.;//Double.parseDouble(Soil.get(INIT_VALUE_SOIL)[1]);//Soil[10][2]
		double Depth_of_soil=Root_Depth*10;//Double.parseDouble(Soil.get(DEPTH_OF_SOIL)[1]);//Soil[11][2]
		String Land_use=Soil.get(LAND_USE)[1];//Soil[12][2]
		double Minimum_LAI=Double.parseDouble(Soil.get(MINIMUM_LAI)[1]);//Soil[13][2]
		double Maximum_LAI=Double.parseDouble(Soil.get(MAXIMUM_LAI)[1]);//Soil[14][2]
		double vcover=Double.parseDouble(Soil.get(INTC_COVERED_FRACTION)[1]);
		double layer_thickness=Double.parseDouble(Soil.get(INTC_LAYER_THICKNESS)[1]);
		double intc_drainage_coeff_b=Double.parseDouble(Soil.get(INTC_DRAINAGE_EXP_B)[1]);
		double intc_drainage_max_d=Double.parseDouble(Soil.get(INTC_DRAINAGE_MAX)[1]);
		double Direct_runoff_factor=Double.parseDouble(Soil.get(DIRECT_RUNOFF_FACTOR)[1]);//Soil[15][2]			
		double Glugla_Coeff_c=Double.parseDouble(Soil.get(GLUGLA_C)[1]);//Soil[16][2]
		double Lambda=Glugla_Coeff_c/Depth_of_soil/Depth_of_soil;///Double.parseDouble(Soil.get(LAMBDA)[1]);//Soil[17][2]
		double Cap_Litter=Double.parseDouble(Soil.get(CAP_LITTER)[1]);//Soil[18][2]
		double Init_Value_Litter=Double.parseDouble(Soil.get(INIT_VALUE_LITTER)[1]);//Soil[19][2]
		double Litter_Reduction_factor=Double.parseDouble(Soil.get(LITTER_REDUCTION_FACTOR)[1]);//Soil[20][2]

		// KN 5/7/23, added a timestep. Date values in Global_Input can be
		// 02.01.1995.00, 02.01.1995.12 
		// ETP_COEFF requires the month to be the second item in the date (split by '.')
		double timestep=Double.parseDouble(Soil.get(TIMESTEP)[1]);  
		
		int landuseIndex = landusesIndexes.get(Land_use); 



		// Check parameters (just to show what is mandatory and what is derived, can be deleted)
		System.out.println("SIMPEL");
		System.out.println("======");
		System.out.println("\nParameters:");
		System.out.println("Time step: "+String.format("%,.2f", timestep)+" hours");
		System.out.println("Field Capacity: "+String.format("%,.2f", Field_Capacity_Percent)+"%");
		System.out.println("Permanent Wilting Point: "+String.format("%,.2f", Permanent_Wilting_Point)+"%");
		System.out.println("Start of Reduction: "+String.format("%,.2f", Start_of_Reduction_Percent)+"%");
		System.out.println("Root depth: "+String.format("%,.2f", Root_Depth)+" cm");
		System.out.println("Init-Value Soil: "+String.format("%,.2f", Init_Value_Soil_Percent)+"%");
		System.out.println("Land use: "+Land_use);
		System.out.println("Minimum_LAI: "+String.format("%,.2f", Minimum_LAI)+"");
		System.out.println("Maximum_LAI: "+String.format("%,.2f", Maximum_LAI)+"");
		System.out.println("Vegetation cover fraction: "+String.format("%,.2f", vcover)+"");
		System.out.println("Layer Thickness: "+String.format("%,.2f", layer_thickness)+"mm");
		System.out.println("Drainage Coeff. b: "+String.format("%,.2f", intc_drainage_coeff_b)+"");
		System.out.println("Max. Drainage rate: "+String.format("%,.2f", intc_drainage_max_d)+"mm/d");
		System.out.println("Capacity of Litter Layer: "+String.format("%,.2f", Cap_Litter)+" mm");
		System.out.println("Initial value of Litter Layer: "+String.format("%,.2f", Init_Value_Litter)+" mm");
		System.out.println("Litter Reduction factor "+String.format("%,.2f", Litter_Reduction_factor)+"");
		System.out.println("Direct runoff factor: "+String.format("%,.2f", Direct_runoff_factor)+"");
		System.out.println("Coefficient c (Glugla): "+String.format("%,.2f", Glugla_Coeff_c)+"");

		System.out.println("\nDerived values");
		System.out.println("Field Capacity: "+String.format("%,.2f", Field_Capacity)+" mm");
		System.out.println("Perm_Wilting_Point: "+String.format("%,.2f", Perm_Wilting_Point)+" mm");
		System.out.println("FWC: "+String.format("%,.2f", fwc)+" mm");
		System.out.println("Start_of_Reduction: "+String.format("%,.2f", Start_of_Reduction)+" mm");
		System.out.println("Init-Value Soil: "+String.format("%,.2f", Init_Value_Soil)+" mm");
		System.out.println("Depth of soil: "+String.format("%,.2f", Depth_of_soil)+" mm");
		System.out.println("Lambda (Glugla): "+String.format("%,.8f", Lambda)+"");

//		# SIMPLE_function
//		# original Excel spreadsheet developed by Georg Hörmann
//		# extended (snow, surface runoff) by Kristian Förster 2022
//		# translated to R by Zoe Bovermann 2023
//      # translated to Java by Kerry Nice 20 June 2023

//		  # prepare bucket model matrix
		double[][] bucket_model = new double[Input.size()][31];
		bucket_model_dates = new ArrayList<String>();
		
//		  # copy Input into bucket-model
//		  #col 1: A Date
//		  # col 2: B Precipitation	
		  for (int i=0;i<Input.size();i++)
		  {
			  bucket_model_dates.add(Input.get(i)[INPUT_DATE]);
			  bucket_model[i][PRECIPITATION]=Double.parseDouble(Input.get(i)[INPUT_P]);
			  if (hasKdown)
			  {
				  bucket_model[i][K_DOWN]=Double.parseDouble(Input.get(i)[INPUT_K_DOWN]);
			  }
			  if (hasWindSpeed)
			  {
				  bucket_model[i][WIND_SPEED]=Double.parseDouble(Input.get(i)[INPUT_WIND_SPEED]);
			  }
			  if (hasLdown)
			  {
				  bucket_model[i][L_DOWN]=Double.parseDouble(Input.get(i)[INPUT_L_DOWN]);
			  }
		  }
		  
		  
//		  # water balance check
		  double sum_prec   = 0.;
		  double sum_etr    = 0.;
		  double sum_runoff = 0.;
		  double init_stor  = Init_Value_Soil; 
		  double init_swe   = 0.;
		    
		  double[][] laiModel = new double[4][4];
		  for (int i=0;i<LAI_model.size();i++)
		  {
			  laiModel[i][LAI_IND] = Double.parseDouble(LAI_model.get(i)[0]);
			  laiModel[i][LAI_DOY] = Double.parseDouble(LAI_model.get(i)[1]);
			  laiModel[i][LAI_PHASE] = Double.parseDouble(LAI_model.get(i)[2]);
		  }
		  
//		  # update LAI model according to soil physics ($kf 2023-03-28) -> ignore land use table
//		  # update Litter according to land use table
		  laiModel[0][LAI_LAI]= Double.parseDouble(Soil.get(MINIMUM_LAI)[1]);
		  laiModel[1][LAI_LAI]= Double.parseDouble(Soil.get(MAXIMUM_LAI)[1]);
		  laiModel[2][LAI_LAI]= Double.parseDouble(Soil.get(MAXIMUM_LAI)[1]);
		  laiModel[3][LAI_LAI]= Double.parseDouble(Soil.get(MINIMUM_LAI)[1]);
//		  # update Litter according to land use table
		  double landuse_DayDegree = Double.parseDouble(Landuse.get(16)[landuseIndex]);
		  
		  double laimodel13 = laiModel[0][LAI_LAI]; //Double.parseDouble(LAI_model.get(0)[2]); // LAI_model[1,3]
		  double laimodel12 = laiModel[0][LAI_DOY]; //Double.parseDouble(LAI_model.get(0)[1]); // LAI_model[1,2]
		  
		  double laimodel23 = laiModel[1][LAI_LAI]; //Double.parseDouble(LAI_model.get(1)[2]); // LAI_model[2,3]
		  double laimodel22 = laiModel[1][LAI_DOY]; //Double.parseDouble(LAI_model.get(1)[1]); // LAI_model[2,2]
		  
		  double laimodel33 = laiModel[2][LAI_LAI]; //Double.parseDouble(LAI_model.get(2)[2]); // LAI_model[3,3]
		  double laimodel32 = laiModel[2][LAI_DOY]; //Double.parseDouble(LAI_model.get(2)[1]); // LAI_model[3,2]
		  
		  double laimodel43 = laiModel[3][LAI_LAI]; //Double.parseDouble(LAI_model.get(3)[2]); // LAI_model[4,3]
		  double laimodel42 = laiModel[3][LAI_DOY]; //Double.parseDouble(LAI_model.get(3)[1]); // LAI_model[4,2]

		  double ts = timestep * 3600; // timestep from soil physics is given in hours, ts is a conversion to seconds
          	  double nt = ts / SECONDS_PER_DAY; // fractional time step
		  // In the model the drying is controlled by the ”litter reduction factor” which 
		  // specifies the maximum evaporation expressed as part of the water content. 
		  // Therefore, a factor 2 indicates that within each step of calculation a maximum 
		  // of half of the storage capacity can evaporate. (source: SIMPEL dcos)
		  // The input value refers to daily time step. Thus it's adjusted to arbitrary
		  // time steps through division by nt.
		  double Litter_Reduction_factor_dt = Litter_Reduction_factor / nt;
		  
		   
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

		    // NEW INTERCEPTION MODEL (which works with different time steps)
		    // Rutter, A.J., Kershaw, K.A., Robins, P.C., Morton, A.J., 1971. 
		    // A predictive model of rainfall interception in forests, 1. Derivation 
		    // of the model from observations in a plantation of Corsican pine. 
		    // Agr. Meteorol. 9, 367–384.

//		    # col 8: H I-Cap
		    bucket_model[krow][I_CAP] = layer_thickness*bucket_model[krow][LAI];

		    double intc_drainage_max = intc_drainage_max_d / 86400 * ts; // max. drainage per time step
		    // Maximum Drainage rate adjusted to capacity
	 	    if(intc_drainage_max > bucket_model[krow][I_CAP]) intc_drainage_max = bucket_model[krow][I_CAP];
		    
		    // define empty storage for first time step and correct for adjustments in C
		    double t_surplus = 0.;
                    double ci = 0.; // first time step: initialize interception capacity, assuming empty storage
		     if(krow>0){
		    	ci = bucket_model[krow-1][I_BAL]; // intercepted water from previous time step
			// Adjust storage as a result of changing capacity
		    	if( bucket_model[krow-1][I_CAP] !=  bucket_model[krow][I_CAP]  && ci >  bucket_model[krow][I_CAP] ){
			        t_surplus = ci - bucket_model[krow][I_CAP];
        			ci = bucket_model[krow][I_CAP];
		        }
    		    }
		    //double ntf = 0.25; // fraction of precipitation that always becomes throughfall, todo: make parameter adjustable
		    double direct_throughfall =  bucket_model[krow][SNOW_MELT_RAIN] *(1.-vcover);
		    double intc_in =  bucket_model[krow][SNOW_MELT_RAIN] - direct_throughfall;
		    double intc_stor_guess = intc_in+ci; // storage depth

//		    # col 9: I ETi: Evaporation from wetted leaf 
		    // bucket_model[krow][INT_ETI_LEAF] = Math.min(bucket_model[krow][I_CAP], bucket_model[krow][ETP_INPUT]);
		    if(intc_stor_guess > bucket_model[krow][I_CAP])
			    bucket_model[krow][INT_ETI_LEAF] = Math.min(bucket_model[krow][ETP_INPUT],intc_stor_guess); // real evaporation
    		    else
			    bucket_model[krow][INT_ETI_LEAF] = Math.min(bucket_model[krow][ETP_INPUT]*intc_stor_guess/ bucket_model[krow][I_CAP], intc_stor_guess);
//		    # col 10: J I-Bal.: Temp calcalation
		    //bucket_model[krow][I_BAL] = bucket_model[krow][SNOW_MELT_RAIN] - bucket_model[krow][INT_ETI_LEAF];
		    double intc_drainage;
		    if(bucket_model[krow][I_CAP] < intc_stor_guess-bucket_model[krow][INT_ETI_LEAF]) 
			intc_drainage = Math.max(intc_stor_guess-bucket_model[krow][INT_ETI_LEAF]-bucket_model[krow][I_CAP],intc_drainage_max);
    		    else
		    	intc_drainage = Math.min(intc_drainage_max*Math.exp(intc_drainage_coeff_b*(intc_stor_guess-bucket_model[krow][INT_ETI_LEAF] - bucket_model[krow][I_CAP])/bucket_model[krow][I_CAP]), intc_stor_guess-bucket_model[krow][INT_ETI_LEAF]);

		    // updated meaning of col 10 "I-BAL" is interception storage
		    bucket_model[krow][I_BAL] = intc_stor_guess-bucket_model[krow][INT_ETI_LEAF]-intc_drainage;
		    if(bucket_model[krow][I_BAL] <0){
      			intc_drainage = Math.max(intc_stor_guess+ bucket_model[krow][I_BAL],0);
      			 bucket_model[krow][I_BAL] = 0.;
    		    }
		    
//		    # col 11: K I-Prec.: Throughfall   
		    bucket_model[krow][I_PREC] = direct_throughfall + intc_drainage + t_surplus;
		    
//		    # col 12: L I-Rem.: Remaining ETa passed to subsequent model
		    // bucket_model[krow][I_REM] = -1*Math.min(bucket_model[krow][I_BAL],0)
		    // 		-bucket_model[krow][INT_ETI_LEAF]+bucket_model[krow][ETP_INPUT];
		    bucket_model[krow][I_REM] = bucket_model[krow][ETP_INPUT] - bucket_model[krow][INT_ETI_LEAF];	

//		    # col 13+14+15
//		    # first row
//		    # litter reduction is now adjusted to time step of input series
		    if(krow==0)
		    {
//		      # col 13: M ETi Litter
		      double[] temp = new double[]{Cap_Litter,(0+bucket_model[krow][I_PREC])/Litter_Reduction_factor_dt};
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
		      double[] temp = new double[]{Cap_Litter, (bucket_model[krow-1][CONTENT]+bucket_model[krow][I_PREC])/Litter_Reduction_factor_dt};
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
//		    # Groundwater recharge computation with Glugla approach now adjusted to time step length
		    if(krow==0)
		    {
//		      # col 17: Q Inf-Limit
		      bucket_model[krow][INF_LIMIT] = (Field_Capacity-Init_Value_Soil)*0.25; // *(1.-Direct_runoff_factor/100.); /* bug fix */
		      
//		      # col 18: R P-Inf
		      bucket_model[krow][P_LINF] = Math.min(bucket_model[krow][S_RESTN],bucket_model[krow][INF_LIMIT])*(1.-Direct_runoff_factor/100.) * nt;
		      
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
		        bucket_model[krow][SEEPAGE] =  Lambda*Math.pow( (bucket_model[krow][ET_BALANCE]-Perm_Wilting_Point),2) * nt;
		      }
		      else
		      {
		        bucket_model[krow][SEEPAGE] =  Lambda*Math.pow( (Field_Capacity-Perm_Wilting_Point),2) * nt;
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
		        		Math.pow((bucket_model[krow][ET_BALANCE]-Perm_Wilting_Point),2) * nt;
		      }
		      else
		      {
		        bucket_model[krow][SEEPAGE] =  Lambda*Math.pow((Field_Capacity-Perm_Wilting_Point),2) * nt;
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
		  double water_balance = sum_prec - sum_etr - sum_runoff + init_swe + init_stor - bucket_model[krowEnd][STORAGE] - bucket_model[krowEnd][SNOW_WATER_EQUI] - bucket_model[krowEnd][I_BAL] - bucket_model[krowEnd][CONTENT];
		  System.out.println("water balance check");
		  System.out.println(water_balance);
		  
		  return bucket_model;
		}
	

	//match input headers to forcing fields and figure out which are in the input. Some are mandatory, some optional
	public void matchInputHeaders(ArrayList<String[]> InputHeader)
	{
		for (int i=0;i<InputHeader.size();i++)
		{
			if (InputHeader.get(i).equals("Date"))
			{
				INPUT_DATE=i;
			}
			else if (InputHeader.get(i).equals("DOY"))
			{
				INPUT_DOY=i;
			}
			else if (InputHeader.get(i).equals("P"))
			{
				INPUT_P=i;
			}
			else if (InputHeader.get(i).equals("T14"))
			{
				INPUT_T14=i;
			}
			else if (InputHeader.get(i).equals("RF14"))
			{
				INPUT_R14=i;
			}
			else if (InputHeader.get(i).equals("KDOWN"))
			{
				INPUT_K_DOWN=i;
				hasKdown = true;
			}
			else if (InputHeader.get(i).equals("LDOWN"))
			{
				INPUT_L_DOWN=i;
				hasLdown = true;
			}
			else if (InputHeader.get(i).equals("WindSpeed"))
			{
				INPUT_WIND_SPEED=i;
				hasWindSpeed = true;
			}
			else if (InputHeader.get(i).equals("ET0"))
			{
				INPUT_ET0=i;
				readETfromFile = true;
			}
		}
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
   public ArrayList<String[]> readLargerTextFileAlternateToArray(String aFileName, boolean skipHeader, String deliminator, boolean onlyFirstLine) 
   {
	   ArrayList<String[]> returnFile = new ArrayList<String[]>();
	   int count = 0;
	    Path path = Paths.get(aFileName);
	    try (BufferedReader reader = Files.newBufferedReader(path, ENCODING))
	    {
	      String line = null;
	      while ((line = reader.readLine()) != null) 
	      {
	    	  if (onlyFirstLine && count > 0)
	    	  {
	    		  return returnFile;
	    	  }
	    	  
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
   
   public TreeMap<String,String[]> readLargerTextFileAlternateToTreeMap(String aFileName, boolean skipHeader, String deliminator) 
   {
	   TreeMap<String,String[]> returnFile = new TreeMap<String,String[]>();
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
	    	  String key = lineSplit[0];
	    	  returnFile.put(key,lineSplit);
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
				+ "ETaTotal" + tab
				+ "K_DOWN" + tab
				+ "WINDSPEED" + tab
				+ "L_DOWN";
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
