package thud.simpel;

public class MaespaCom
{
	
	
    public static int KHRS = 48;                  // Number of time intervals in a day
    public static double HHRS = (KHRS) / 2.0;           // Half a day length
    public static double SPERHR = 3600.0 * 24.0 / KHRS; // Seconds in one time interval

    // Values of physical constants
    public final static double DEFWIND = 2.5;       // Default wind speed (m s-1)
    public final static double UMOLPERJ = 4.57;     // Conversion from J to umol quanta
    public final static double FPAR = 0.5;          // Fraction of global radiation that is PAR
    public final static double ABSZERO = -273.15;   // Absolute zero in degrees Celsius
    public final static double FREEZE = 273.15;     // Zero degrees Celsius in Kelvin
    public final static double TAU = 0.76;          // Transmissivity of atmosphere
    public final static double PI = 3.1415927;      // Pi
    public final static double TWOPI = 2.0 * PI;    // Two times Pi
    public final static double PID2 = PI / 2.0;     // Pi divided by two
    public final static double PID180 = PI / 180.0; // Pi divided by 180 degrees
    public final static double AIRMA = 29.e-3;      // mol mass air (kg/mol)
    public final static double PATM = 1.0125E5;     // atmospheric pressure - standard condns (Pa)
    public final static double CPAIR = 1010.0;      // heat capacity of air (J kg-1 K-1)
    public final static double CPH2O = 4.186E06;    // heat capacity of water (J kg-1 K-1)
    public final static double CPQUARTZ = 1.942E06; // heat capacity of quartz (J kg-1 K-1)
    public final static double TCQUARTZ = 7.7;      // thermal conductivity of quartz (W m-1 K-1)
    public final static double TCH2O = 0.594;       // thermal conductivity of water (W m-1 K-1)
    public final static double TCORG = 0.25;        // thermal conductivity of organic matter (W m-1 K-1)
    public final static double SOILALBEDO = 0.15;   // Albedo of soil, without snow.
    public final static double DHEAT = 21.5e-6;     // molecular diffusivity for heat
    public final static double EMLEAF = 0.95;       // Emissivity of thermal radiation by leaf
    public final static double EMSOIL = 0.95;       // Emissivity of thermal radiation by soil
    public final static double H2OLV0 = 2.501e6;    // latent heat H2O (J/kg)
    public final static double H2OMW = 18.e-3;      // mol mass H2O (kg/mol)
    public final static double H2OVW = 18.05e-6;    // partial molal volume of water at 20C (m3 mol-1)
    public final static double RCONST = 8.314;      // universal gas constant (J/mol/K)
    public final static double SIGMA = 5.67e-8;     // Steffan Boltzman constant (W/m2/K4)
    public final static double GBHGBC = 1.32;       // Ratio of Gbh:Gbc
    public final static double GSVGSC = 1.57;       // Ratio of Gsw:Gsc
    public final static double GBVGBH = 1.075;      // Ratio of Gbw:Gbh
    public final static double ALPHAQ = 0.425;      // Quantum yield of RuBP regen (mol mol-1)
    public final static double SOLARC = 1370;       // Solar constant (J m-2 s-1)
    public final static double GCPERMOL = 12.0;     // Grams // per mol //
    public final static double CPERDW = 0.5;        // fraction per DW
    public final static double VONKARMAN = 0.41;    // von Karman's constant
    public final static double GRAV = 9.8067;       // Gravitational acceleration
    
}
