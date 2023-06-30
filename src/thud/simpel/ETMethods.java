package thud.simpel;


public class ETMethods
{

	public static void main(String[] args)
	{
		

	}
	
	//unit tests to demonstrate how the methods should be called
	public void testetpot()
	{		
		double laiday=0.000000E+00;
		double albday=1.000000E-01;
		double cht=0.000000E+00;
		double co2=3.300000E+02;
		double gsi=0.000000E+00;
		double hru_ra=4.259289E+00;
		double hru_rmx=1.276488E+01;
		double hru_sub=1;
		double icr=1;
		double idplt=98;
		double igro=0;
		int ihru=1;
		int ipet=0;
		double nro=1;
		double petmeas=0.000000E+00;
		double rhd=9.900000E-01;
		double sno_hru=0.000000E+00;
		double sub_elev=3.990000E+02;
		double tmn=2.000000E+00;
		double tmpav=1.900000E+01;
		double tmx=3.600000E+01;
		double u10=0.000000E+00;
		double vpd2=0.000000E+00;
		double harg_petco=2.300000E-03;
		double[] returnValue = etpot(laiday, albday, cht, co2, gsi, hru_ra, hru_rmx, hru_sub, icr, idplt, igro, 
				ihru, ipet, nro, petmeas, rhd, sno_hru, sub_elev, tmn, tmpav, tmx, u10, vpd2, harg_petco);
		
		double[] etpotExpected =  {0.000000E+00,        4.798565E-01,        2.198601E-02};		
//	    assertArrayEquals(etpotExpected, returnValue, 0.001);
	    
		laiday=0.000000E+00;
		albday=1.000000E-01;
		cht=0.000000E+00;
		co2=3.300000E+02;
		gsi=0.000000E+00;
		hru_ra=2.593315E+00;
		hru_rmx=1.285207E+01;
		hru_sub=1;
		icr=1;
		idplt=40;
		igro=1;
		ihru=14;
		ipet=0;
		nro=1;
		petmeas=0.000000E+00;
		rhd=9.870892E-01;
		sno_hru=0.000000E+00;
		sub_elev=3.990000E+02;
		tmn=2.000000E+00;
		tmpav=1.900000E+01;
		tmx=3.600000E+01;
		u10=0.000000E+00;
		vpd2=0.000000E+00;
		harg_petco=2.300000E-03;
		returnValue = etpot(laiday, albday, cht, co2, gsi, hru_ra, hru_rmx, hru_sub, icr, idplt, igro, 
				ihru, ipet, nro, petmeas, rhd, sno_hru, sub_elev, tmn, tmpav, tmx, u10, vpd2, harg_petco);		
		double[] etpotExpected2 =  {0.000000E+00,
		        2.263287E-01,
		        2.838588E-02};		
//	    assertArrayEquals(etpotExpected2, returnValue, 0.001);   		    
	}
	
	
	public void testWETEVAP()
	{
		double WIND=2.9000001;
		double ZHT=30;
		double Z0HT=3;
		double ZPD=16;
		double PRESSPA=88100;
		double TAIRC=12.1400146;
		double RNET=0;
		double VPDPA=1308.64392;
		double CANSTORE=0;
		double MAXSTORAGE=0.5;
		double EVAPSTORE=Double.NaN;
		double PPT=0;
		double WETEVAP = WETEVAP(WIND, ZHT, Z0HT, ZPD, PRESSPA, TAIRC, RNET, VPDPA, CANSTORE, MAXSTORAGE, EVAPSTORE, PPT);
		double WETEVAPExpected =0 ;
//		assertEquals(WETEVAPExpected, WETEVAP, 0.001);
	}
	
	public void testETCAN()
	{
		double WIND=2.9000001;
		double ZHT=30;
		double Z0HT=3;
		double ZPD=16;
		double PRESS=88100;
		double TAIR=12.1400003;
		double RNET=0;
		double VPD=1308.64392;
		double GSCAN=0;
		double STOCKING=0.000600000028;
		double ETCAN = ETCAN(WIND, ZHT, Z0HT, ZPD, PRESS, TAIR, RNET, VPD, GSCAN, STOCKING);
		double ETCANExpected =0 ;
//		assertEquals(ETCANExpected, ETCAN, 0.001);				               

		 WIND=4.770000;
		 ZHT=4.000000 ;
		 Z0HT=2.0000000E-02;
		 ZPD=6.6000000E-02;
		 PRESS=996.3600 ;
		 TAIR=11.95999 ;
		 RNET=0.0000000E+00;
		 VPD=362.1992;
		 GSCAN= 1.0000000E+09;
		 STOCKING=1.000000  ;
		 ETCAN = ETCAN(WIND, ZHT, Z0HT, ZPD, PRESS, TAIR, RNET, VPD, GSCAN, STOCKING);
		 ETCANExpected =30.80120  ;
//		assertEquals(ETCANExpected, ETCAN, 0.4);		
		
		 WIND=2.9000001;
		 ZHT=30;
		 Z0HT=3;
		 ZPD=16;
		 PRESS=88100;
		 TAIR=12.1400146;
		 RNET=0;
		 VPD=1308.64392;
		 GSCAN=1e+09;
		 STOCKING=1;
		 ETCAN = ETCAN(WIND, ZHT, Z0HT, ZPD, PRESS, TAIR, RNET, VPD, GSCAN, STOCKING);
		 ETCANExpected =43297.4062 ;
//		assertEquals(ETCANExpected, ETCAN, 0.4);
	}
	
	public void testGBCAN()
	{
		double WIND=2.9000001;
		double ZHT=30;
		double Z0HT=3;
		double ZPD=16;
		double PRESS=88100;
		double TAIR=12.1400003;
		double GBCAN = GBCAN(WIND, ZHT, Z0HT, ZPD, PRESS, TAIR);
		double GBCANExpected =7.63049412 ;
//		assertEquals(GBCANExpected, GBCAN, 0.001);
		
		 WIND=2.9000001;
		 ZHT=30;
		 Z0HT=3;
		 ZPD=16;
		 PRESS=88100;
		 TAIR=12.1400146;
		 GBCAN = GBCAN(WIND, ZHT, Z0HT, ZPD, PRESS, TAIR);
		 GBCANExpected =7.63049412 ;
//		assertEquals(GBCANExpected, GBCAN, 0.001);	
	}
	
	public void testPENMON()
	{
		double PRESS=88100;
		double SLOPE=93.8098145;
		double LHV=44501.1953;
		double RNET=0;
		double VPD=1308.64392;
		double GH=7.63049412;
		double GV=7.63049459;
		double PENMON = PENMON(PRESS, SLOPE, LHV, RNET, VPD, GH, GV);
		double PENMONExpected =0.0432974063;
//		assertEquals(PENMONExpected, PENMON, 0.001);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	//method translated from FORTRAN SWAT model
	/** this subroutine calculates potential evapotranspiration using one of three methods. If Penman-Monteith is being used, potential plant transpiration is also calculated.
	 * @param laiday leaf area index (m**2/m**2 )
	 * @param albday albedo for the day in HRU
	 * @param cht canopy height (m)
	 * @param co2 CO2 concentration (ppmv)
	 * @param gsi maximum stomatal conductance (m/s)
	 * @param hru_ra solar radiation for the day in HRU (MJ/m^2)
	 * @param hru_rmx maximum possible radiation for the day in HRU (MJ/m^2)
	 * @param hru_sub subbasin in which HRU is located
	 * @param icr sequence number of crop grown within the current year
	 * @param idplt land cover code from crop.dat
	 * @param igro land cover status code (0 no land cover currently growing, 1 land cover growing)
	 * @param ihru HRU numbe
	 * @param ipet code for potential ET method (0 Priestley-Taylor method, 1 Penman/Monteith method, 2 Hargreaves method, 3 read in PET values)
	 * @param nro sequence number of year in rotation
	 * @param petmeas potential ET value read in for day (mm H2O)
	 * @param rhd relative humidity for the day in HRU
	 * @param sno_hru amount of water in snow in HRU on current day (mm H2O)
	 * @param sub_elev elevation of HRU (m)
	 * @param tmn minimum air temperature on current day in HRU (deg C)
	 * @param tmpav average air temperature on current day in HRU (deg C)
	 * @param tmx maximum air temperature on current day for HRU (deg C)
	 * @param u10 wind speed (measured at 10 meters above surface (m/s) 
	 * @param vpd2 rate of decline in stomatal conductance per unit increase in vapor pressure deficit ( (m/s)*(1/kPa) )
	 * @param harg_petco coefficient related to radiation used in hargreaves eq (range: 0.0019 - 0.0032)
	 * @return
	 */
	public double[] etpot(double laiday, double albday, double cht, double co2, double gsi, double hru_ra, double hru_rmx, double hru_sub,
			          double icr, double idplt, double igro, int ihru, int ipet, double nro, double petmeas, double rhd,
			          double sno_hru, double sub_elev, double tmn, double tmpav, double tmx, double u10, double vpd2, double harg_petco)
	{	
		
		//return values
		double ep_max =	0 ; //maximum amount of transpiration (plant et) that can occur on current day in HRU (mm H2O)
		double pet_day =0; // potential evapotranspiration on current day in HRU (mm H2O)
		double vpd; //vapor pressure deficit (kPa)
		//end return values
		
		//int  j; // HRU number
		double tk; //average air temperature on current day for HRU
		double pb; //mean atmospheric pressure (kPa)
		double gma; //psychrometric constant (kPa/deg C)
		double xl; //latent heat of vaporization (MJ/kg)
		double ea; //saturated vapor pressure (kPa)
		double ed; //actual vapor pressure (kPa)
		double dlt; //slope of the saturation vapor pressure-temperature curve (kPa/deg C)
		double ramm; //extraterrestrial radiation (MJ/m2)
		double ralb1; //net incoming radiation (MJ/m2)
		double ralb; //net incoming radiation for PET (MJ/m2)
		double xx; //difference between vpd and vpthreshold (kPa)
		double rbo; //net emissivity
		double rto; //cloud cover factor
		double rn; //net radiation (MJ/m2)
		double uzz; //wind speed at height zz (m/s)
		double zz; //height at which wind speed is determined (cm)
		double zom; //roughness length for momentum transfer (cm)
		double zov; //roughness length for vapor transfer (cm)
		double rv; //aerodynamic resistance to sensible heat and vapor transfer (s/m)
		double rn_pet; //net radiation for continuous crop cover (MJ/m2)
		double fvpd; //amount of vapro pressure deficit over threshold value (kPa)
		double rc; //canopy resistance (s/m)
		double rho; //K1*0.622*xl*rho/pb ( MJ/(m3*kPa) )
		double rout; //outgoing radiation (MJ/m2)
		double d; //displacement height for plant type (cm)
		double chz; //vegetation height (cm)
		double gsi_adj;
		double pet_alpha; //alpha factor in Priestley-Taylor PET equation

		tk = 0.;
		tk = tmpav + 273.15;

		// calculate mean barometric pressure
		pb = 0.;
		pb = 101.3 - sub_elev *  (0.01152 - 0.544e-6 * sub_elev);

		// calculate latent heat of vaporization
		xl = 0.;
		xl = 2.501 - 2.361e-3 * tmpav;

		// calculate psychrometric constant
		gma = 0.;
		gma = 1.013e-3 * pb / (0.622 * xl);

		// calculate saturation vapor pressure, actual vapor pressure and
		// vapor pressure deficit
		ea = 0.;
		ed = 0.;
		vpd = 0.;
		ea = ee(tmpav);
		ed = ea * rhd;
		vpd = ea - ed;

		//calculate the slope of the saturation vapor pressure curve
		dlt = 0.;
		dlt = 4098. * ea / Math.pow((tmpav + 237.3), 2);
		    
		// DETERMINE POTENTIAL ET
		switch (ipet)
		{
			case 0:  // PRIESTLEY-TAYLOR POTENTIAL EVAPOTRANSPIRATION METHOD	
				// net radiation
		        // calculate net short-wave radiation for PET
		        ralb = 0.;
		        if (sno_hru <= .5) 
		        {
		            ralb = hru_ra * (1.0 - 0.23);
		        }
		        else
		        {
		        	ralb = hru_ra * (1.0 - 0.8);
		        }
		        // calculate net long-wave radiation
		        // net emissivity  equation 2.2.20 in SWAT manual
		        rbo = 0.;
		        rbo = -(0.34 - 0.139 * Math.sqrt(ed) );

		        // cloud cover factor equation 2.2.19
		        rto = 0.;
		        if (hru_rmx < 1.e-4) 
		        {
		            rto = 0.;
		        }
		        else
		        {
		            rto = 0.9 * (hru_ra / hru_rmx) + 0.1;
		        }
		
	            // net long-wave radiation equation 2.2.21
		        rout = 0.;
		        rout = rbo * rto * 4.9e-9 *  Math.pow(tk, 4) ;

		        // calculate net radiation
		        rn_pet = 0.;
		        rn_pet = ralb + rout;
		        // net radiation

		        pet_alpha = 1.28;
		        pet_day = pet_alpha * (dlt / (dlt + gma)) * rn_pet / xl;
		        pet_day = Math.max(0., pet_day);
								
				break;
				
			case 1: // PENMAN-MONTEITH POTENTIAL EVAPOTRANSPIRATION METHOD
				// net radiation
			    // calculate net short-wave radiation for PET
			    ralb = 0.;
			    if (sno_hru <= .5) 
			    {
			    	ralb = hru_ra * (1.0 - 0.23) ;
			    }
			    else
			    {
			        ralb = hru_ra * (1.0 - 0.8) ;
			    }
			         
			    // calculate net short-wave radiation for max plant ET
			    ralb1 = 0.;
			    ralb1 = hru_ra * (1.0 - albday) ;

			    // calculate net long-wave radiation
			    // net emissivity  equation 2.2.20 in SWAT manual
			    rbo = 0.;
			    rbo = -(0.34 - 0.139 * Math.sqrt(ed));

			    // cloud cover factor equation 2.2.19
			    rto = 0.;
			    if (hru_rmx < 1.e-4) 
			    {
			    	rto = 0.;
			    }			            
			    else
			    {
			    	rto = 0.9 * (hru_ra / hru_rmx) + 0.1;
			    }
			              
			    // net long-wave radiation equation 2.2.21
			    rout = 0.;
			    rout = rbo * rto * 4.9e-9 * Math.pow(tk, 4) ;

			    // calculate net radiation
			    rn = 0.;
			    rn_pet = 0.;
			    rn = ralb1 + rout;
			    rn_pet = ralb + rout;
			    // net radiation

			    rho = 0.;
			    rho = 1710. - 6.85 * tmpav;

			    if (u10 < 0.01) 
			    {
			    	u10 = 0.01;
			    }

			    // potential ET: reference crop alfalfa at 40 cm height
			    rv = 0.;
			    rc = 0.;
			    rv = 114. / (u10 * Math.pow((170./1000.), 0.2) );
			    rc = 49. / (1.4 - 0.4 * co2 / 330.);
			    pet_day = (dlt * rn_pet + gma * rho * vpd / rv) / (xl * (dlt + gma * (1. + rc / rv)));

			    pet_day = Math.max(0., pet_day);
			 
			    // maximum plant ET
			    if (igro <= 0) 
			    {
			        ep_max = 0.0;
			    }
			    else
			    {
			    	// determine wind speed and height of wind speed measurement
			        // adjust to 100 cm (1 m) above canopy if necessary
			        uzz = 0.;
			        zz = 0.;
			        if (cht <= 1.0) 
			        {
			              zz = 170.0;
			        }
			        else
			        {
			              zz = cht * 100. + 100.;
			        }
			         
			        uzz = u10 * Math.pow((zz/1000.), 0.2);

			        // calculate canopy height in cm
			        chz = 0.;
			        if (cht < 0.01) 
			        {
			              chz = 1.;
			        }
			        else
			        {
			              chz = cht * 100.;
			        }
			       
			        // calculate roughness length for momentum transfer
			        zom = 0.;
			        if (chz <= 200.) 
			        {
			              zom = 0.123 * chz;
			        }
			        else
			        {
			              zom = 0.058 * Math.pow(chz, 1.19) ;
			        }
			      
			        // calculate roughness length for vapor transfer
			        zov = 0.;
			        zov = 0.1 * zom;

			        // calculate zero-plane displacement of wind profile
			        d = 0.;
			        d = 0.667 * chz;

			        // calculate aerodynamic resistance
			        rv = Math.log((zz - d) / zom) * Math.log((zz - d) / zov);
			        rv = rv / ( Math.pow((0.41), 2) * uzz);

			        // adjust stomatal conductivity for low vapor pressure
			        // this adjustment will lower maximum plant ET for plants
			        // sensitive to very low vapor pressure
			        xx = 0.;
			        fvpd = 0.;
			        xx = vpd - 1.;
			        if (xx > 0.0) 
			        {
			              fvpd = Math.max(0.1,1.0 - vpd2 * xx);
			        }
			        else
			        {
			              fvpd = 1.0;
			        }
			   
			        gsi_adj = 0.;
			        gsi_adj = gsi * fvpd;
			            
			        // calculate canopy resistance
			        rc = 0.;
			        rc = 1. / gsi_adj;               //single leaf resistance
			        rc = rc / (0.5 * (laiday + 0.01)  * (1.4 - 0.4 * co2 / 330.));

			        // calculate maximum plant ET
			        ep_max = (dlt * rn + gma * rho * vpd / rv) / (xl * (dlt + gma * (1. + rc / rv)));
			        if (ep_max < 0.) 
			        {
			        	ep_max = 0.;
			        }
			        ep_max = Math.min(ep_max, pet_day);
			          
			    }
				
				break;
			case 2: // HARGREAVES POTENTIAL EVAPOTRANSPIRATION METHOD
				// extraterrestrial radiation
		        // 37.59 is coefficient in equation 2.2.6 !!extraterrestrial
		        // 30.00 is coefficient in equation 2.2.7 !!max at surface
		        ramm = 0.;
		        ramm = hru_rmx * 37.59 / 30. ;

		        if (tmx > tmn) 
		        {
		         pet_day = harg_petco * (ramm / xl)*(tmpav + 17.8)*  Math.pow((tmx - tmn), 0.5);
		         pet_day = Math.max(0., pet_day);
		        }
		        else
		        {
		          pet_day = 0.;
		        }
				
				break;
			case 3: // READ IN PET VALUES
				pet_day = petmeas;
				break;
		       

			default:
				break;
		}
		
		double[] returnValues = {ep_max, pet_day, vpd};
		return returnValues;

	}
	
	/**
	 * This function calculates saturation vapor pressure at a given air
	 * temperature.
	 * 
	 * @param tk  mean air temperature (deg C)
	 * @return saturation vapor pressure (kPa)
	 */

	public double ee(double tk)
	{
		double ee = 0.;
		if (tk + 237.3 != 0.)
		{
			ee = (16.78 * tk - 116.9) / (tk + 237.3);
			ee = Math.exp(ee);
		}
		return ee;

	}
	

//method translated from FORTRAN MAESPA model
    /**  Penman-monteith equation without stomatal limitation
(= Penman equation) for wet canopies.
Modified from SPA version; now calls the maestra function ETCAN
with infinite canopy conductance.
RAD June 2008.
     * @param WIND
     * @param ZHT
     * @param Z0HT
     * @param ZPD
     * @param PRESSPA
     * @param TAIRC
     * @param RNET
     * @param VPDPA
     * @param CANSTORE
     * @param MAXSTORAGE
     * @param EVAPSTORE
     * @param PPT
     * @return
     */
    public double WETEVAP(double WIND, double ZHT, double Z0HT, double ZPD, 
    		double PRESSPA, double TAIRC, double RNET, 
    		double VPDPA, double CANSTORE, double MAXSTORAGE, 
    		double EVAPSTORE, double PPT)
    {
	    double STOCK,GCAN,POTEVAPMUMOL,POTENTIALEVAP,RATIO,ETCAN;
	
	    STOCK = 1.;   //! Dummy (to avoid unit conversion in ETCAN).
	    GCAN = 1E09;  //! an arbitrary large number (~Inf).
	
	    //! Set VPD to near zero if it is raining (otherwise get very high wet evaporation rates).
	    if (PPT > 0.0)
	    {
	    	VPDPA = 1;
	    }
	
	    //! Potential evaporation from a wet canopy in mu mol m-2 s-1.
	    POTEVAPMUMOL = ETCAN(WIND,ZHT,Z0HT,ZPD,PRESSPA,TAIRC,RNET,VPDPA,GCAN,STOCK);
	
	    //! Convert to mm timestep-1.
	    POTENTIALEVAP = POTEVAPMUMOL * MaespaCom.SPERHR * 18 * 1E-09;
	
	    //! Modifier to potential ET: canopy storage / maximum storage (Rutter).
	    //! rate of evaporation from storage nb store cannot exceed max storage
	    RATIO = Math.min(1.,CANSTORE/MAXSTORAGE);
	    EVAPSTORE = POTENTIALEVAP * RATIO;
	
	    return EVAPSTORE;
    }
  

	/**  Calculate transpiration by applying Penman-Monteith to whole canopy.
	 Returns umol m-2 s-1.
	 * @param WIND
	 * @param ZHT
	 * @param Z0HT
	 * @param ZPD
	 * @param PRESS
	 * @param TAIR
	 * @param RNET
	 * @param VPD
	 * @param GSCAN
	 * @param STOCKING
	 * @return
	 */
	public double ETCAN(double WIND, double ZHT, double Z0HT, double ZPD, double PRESS, double TAIR, double RNET, double VPD, double GSCAN, double STOCKING)
	{    
	    double LHV,GB,GSV,RNETM2,SLOPE,GH,GV;
	    double ETCAN;
	
	    //! Get boundary layer conductance
	    GB = GBCAN(WIND,ZHT,Z0HT,ZPD,PRESS,TAIR);
	    
	    //! Convert mol CO2/tree/s to mol H2O/m2/s
	    GSV = GSCAN*MaespaCom.GSVGSC*STOCKING;
	    RNETM2 = RNET*STOCKING;
	
	    if (GB*GSV > 0.0) 
	    {
	        //! Latent heat of water vapour at air temperature (J mol-1)
	        LHV = HEATEVAP(TAIR) * MaespaCom.H2OMW;
	
	        //! Const s in Penman-Monteith equation  (Pa K-1)
	        SLOPE = (SATUR(TAIR + 0.1) - SATUR(TAIR)) / 0.1;
	        
	        //! Call Penman-Monteith
	        GH = GB;
	        GV = 1./(1./GSV + 1./GB);
	        ETCAN = PENMON(PRESS,SLOPE,LHV,RNETM2,VPD,GH,GV)*1E6;
	    }
	    else
	    {
	        ETCAN = 0.0;
	    }
	
	    return ETCAN;
	}



	/**  This subroutine calculates evapotranspiration by leaves using the Penman-Monteith equation.
	 Inputs:      PRESS atmospheric pressure, Pa
	            SLOPE slope of VPD/T curve, Pa K-1
	            LHV latent heat of water at air T, J mol-1
	            RNET net radiation, J m-2 s-1
	            VPD vapour pressure deficit of air, Pa
	            GH boundary layer conductance to heat (free & forced & radiative components), mol m-2 s-1
	            GV conductance to water vapour (stomatal & bdry layer components), mol m-2 s-1
	 Result in mol H2O m-2 s-1.
	 * @param PRESS
	 * @param SLOPE
	 * @param LHV
	 * @param RNET
	 * @param VPD
	 * @param GH
	 * @param GV
	 * @return
	 */
	public double PENMON(double PRESS, double SLOPE, double LHV, double RNET, double VPD, double GH, double GV)
	{
	    double ET,GAMMA;
	    double PENMON;
	    
	    GAMMA = MaespaCom.CPAIR*MaespaCom.AIRMA*PRESS/LHV;
	
	    if (GV > 0.0) 
	    {
	        ET = (SLOPE * RNET + VPD * GH * MaespaCom.CPAIR * MaespaCom.AIRMA) / (SLOPE + GAMMA * GH/GV);
	    }
	    else
	    {
	        ET = 0.0;
	    }
	    PENMON = ET / LHV;
	    //!      if (PENMON < 0.0) PENMON = 0.0            //! BM 12/05 Should not be negative	
	      return PENMON;
	}




	/**  Canopy boundary layer conductance (from Jones 1992 p 68)
	 in mol m-2 s-1
	 ZHT =
	 * @param WIND
	 * @param ZHT
	 * @param Z0HT
	 * @param ZPD
	 * @param PRESS
	 * @param TAIR
	 * @return
	 */
	public double GBCAN(double WIND, double ZHT, double Z0HT, double ZPD, double PRESS, double TAIR)
	{
	    double CMOLAR;
	    double GBCAN;
	    
	    if (Z0HT > 0.0) 
	    {
	        //! Formula from Jones 1992 p 68
	        GBCAN = WIND*( Math.pow(MaespaCom.VONKARMAN,2) )/Math.pow((Math.log((ZHT - ZPD)/Z0HT)),2);
	        
	        //! Convert from m s-1 to mol m-2 s-1
	        CMOLAR = PRESS / (MaespaCom.RCONST * TK(TAIR));
	        GBCAN = GBCAN * CMOLAR;
	    }
	    else
	    {
	        GBCAN = 0.0;
	    }
	
	    return GBCAN;
	}


	/** Converts Celsius temperature to Kelvin.
	 * @param TCELSIUS
	 * @return
	 */
	public double TK( double TCELSIUS)
	{
	//! Converts Celsius temperature to Kelvin.	    
	    double TK = TCELSIUS - MaespaCom.ABSZERO;
	    return TK;
	}

	//!**********************************************************************
	/** Calculate saturated water vapour pressure (Pa) at temperature TAC (Celsius) from Jones 1992 p 110 (note error in a - wrong units)
	* @param TAC
	* @return SATUR
	*/
	public double SATUR(double TAC)
	{
	//! Calculate saturated water vapour pressure (Pa) at temperature TAC (Celsius)
	//! from Jones 1992 p 110 (note error in a - wrong units)
	//!**********************************************************************
	//  IMPLICIT NONE
	//  REAL TAC
	  double SATUR = 613.75*Math.exp(17.502*TAC/(240.97+TAC));
	  return SATUR;
	}


	/**
	 * @param TAIR
	 * @return
	 */
	public double HEATEVAP(double TAIR)
	{
	//! Latent heat of vaporisation of water (J kg-1)
	//! RAD, May 2008.
	//! TAIR is in celsius.
	//    ! H2OLV is a constant in MAESTCOM        
	    double HEATEVAP = (MaespaCom.H2OLV0 - 2.365E3 * TAIR);
	        
	//    ! A warning    
	    if (TAIR>100.0)
	    {
	    	System.out.println("You are passing T in K to HEATEVAP");
	    }
	    return HEATEVAP;	        
	}

}
