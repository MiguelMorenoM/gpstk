/*********************************************************************
*  $Id:$
*
*  Test program from July 2011.  Written to test the CNAVClock 
*  module as it was being developed.
*
*********************************************************************/
#include <stdio.h>

#include "CNAVClock.hpp"
#include "EngEphemeris.hpp"
#include "CivilTime.hpp"
#include "CommonTime.hpp"
#include "GPSWeekSecond.hpp"

using namespace std;
using namespace gpstk;

int main( int argc, char * argv[] )
{
      // Set time to Day 153, 2011 (6/2/2011) at noon
   CivilTime g( 2011, 6, 2, 12, 14, 44.0, TimeSystem::GPS );
   CommonTime dt = g.convertToCommonTime();

      // Test data (copied from navdmp output for .....)
      // Generally, we'd load these data from the file
   std::string SysID = "G";
   ObsID obsID( ObsID::otNavMsg, ObsID::cbL2, ObsID::tcC2LM );
   ObsID obsID2( ObsID::otNavMsg, ObsID::cbL5, ObsID::tcIQ5 );
   short PRNID     = 3;
   double Toc      = 388800.0;
   short TOWWeek   = 1638;     // By rules of Clock Correction, this must be week of Toc
   double accuracy = 10.61;

      // Test Data copied from RINEX file	
   short rTOWWeek   = 1638;    // By rules of Clock Correction, this must be week of Toc
   double raccuracy = 10.61;
   double rToc      = 388800.0;
   double raf0      = 7.23189674318E-04;
   double raf1      = 5.11590769747E-12;
   double raf2      = 0.0;
   long TOWMsg_1    = 382500;
   long Top         = 378000;
   short AlertMsg   = 0;
   short URAoc_1    = 4;
   short URAoc1_1   = 1;
   short URAoc2_1   = 2;

      // Set time to Day 156, 2011 (6/5/2011) at 1 am
   CivilTime ct2( 2011, 6, 5, 1, 0, 0.0, TimeSystem::GPS );
   CommonTime dt2 = ct2.convertToCommonTime();

      // Test data (copied from navdmp output for PRN 7 Day 156, 2011 at 00:00:00 Transmit Time)
      // Generally, we'd load these data from the file
   short PRNID2     = 7;
   double Toc2      = 7200.0;
   short TOWWeek2   = 1639;     // By rules of Clock Correction, this must be week of Toc
   double accuracy2 = 10.61;
   double af0_2     = 1.32815912E-05;
   double af1_2     = 1.25055521E-12;
   double af2_2     = 0.0;
   long TOWMsg_2    = 0;
   long Top2        = 601200;
   short AlertMsg2  = 1;
   short URAoc_2    = 1;
   short URAoc1_2   = 2;
   short URAoc2_2   = 3;


      // Set time to Day 156, 2011 (6/5/2011) at midnight
   CivilTime ct3( 2011, 6, 5, 0, 0, 0.0, TimeSystem::GPS );
   CommonTime dt3 = ct2.convertToCommonTime();

      // Test data (copied from navdmp output for PRN 9 Day 155, 2011 at 22:00:00 Transmit Time)
      // Generally, we'd load these data from the file
   short PRNID3     = 9;
   double Toc3      = 0.0;
   short TOWWeek3   = 1638;     // By rules of Clock Corection, this must be week of Toc
   double accuracy3 = 10.61;
   long TOWMsg_3    = 597600;
   long Top3        = 594000;
   double af0_3     = 8.43554735E-05;
   double af1_3     = 2.38742359E-12;
   double af2_3     = 0.0;
   short health3    = 0;
   long TOW3        = 597600; 
   short URAoc_3    = 1;
   short URAoc1_3   = 2;
   short URAoc2_3   = 3;

  
   long subframe1[10] = { 0x22C2663D, 0x1F0E29B8, 0x2664002B, 0x09FCC1B6, 0x0F60EB8A,
                          0x1299CE93, 0x29CD3DB6, 0x0597BB0F, 0x00000B68, 0x17B28E5C };
   long subframe2[10] = { 0x22C2663D, 0x1F0E4A28, 0x05809675, 0x0EBD8AF1, 0x00089344,
                         0x008081F8, 0x1330CC2C, 0x0461E855, 0x034F8045, 0x17BB1E68 };
   long subframe3[10] = { 0x22C2663D, 0x1F0E6BA0, 0x3FE129CD, 0x26E31837, 0x0006C96A,
                          0x35A74DFC, 0x065C8B0F, 0x1E4F400A, 0x3FE8966D, 0x05860C44 };

   cout.precision(11);

      // First test case.  Create an CC object with data available from RINEX file.
   cout << endl << "Test Case 1: Creating CC object with data from RINEX file." << endl;
   cout << "Time = " << g << endl;
   CNAVClock cc1;
   cc1.loadData( SysID, obsID, PRNID, AlertMsg, TOWMsg_1, rTOWWeek, Top,
                 rToc, raccuracy, URAoc_1, URAoc1_1, URAoc2_1, raf0, raf1, raf2 );

   double ClkCorr1 = cc1.svClockBias( dt );
   double ClkDrift1 = cc1.svClockDrift( dt ); 
   cout << "Clock Bias cc1:         " << ClkCorr1 << endl;
   cout << "Clock Drift cc1:        " << ClkDrift1 << endl;
   cout << "Time of Prediction cc1: " << GPSWeekSecond(cc1.getTimeOfPrediction()).printf("%F, %g") << endl;
   cout << "CNAV Accuracy Test:     " << SV_CNAV_ACCURACY_MAX_INDEX[URAoc_1+15] << endl;
   cout << "legacy Accuracy Test:   " << SV_ACCURACY_MAX_INDEX[URAoc_1] << endl;

      // Second test case.  Create an CC object with data available from navdump.
   cout << endl << "Test Case 2: Creating CC object with data from navdump." << endl;
   cout << "Time = " << ct2 << endl;
   CNAVClock cc2;
   cc2.loadData( SysID, obsID2, PRNID2, AlertMsg2, TOWMsg_2, TOWWeek2, Top2,
                 Toc2, accuracy2, URAoc_2, URAoc1_2, URAoc2_2, af0_2, af1_2, af2_2 );

   double ClkCorr2 = cc2.svClockBias( dt2 );
   double ClkDrift2 = cc2.svClockDrift( dt2 ); 
   cout << "Clock Bias cc2:  "        << ClkCorr2 << endl;
   cout << "Clock Drift cc2: "        << ClkDrift2 << endl; 
   cout << "Time of Prediction cc2: " << GPSWeekSecond(cc2.getTimeOfPrediction()).printf("%F, %g") << endl;

      // Third test case.  Create an CC object with data available from navdump.
   cout << endl << "Test Case 3: Creating CC object with data from navdump." << endl;
   cout << "Time = " << ct3 << endl;
   CNAVClock cc3;
   cc3.loadData( SysID, obsID, PRNID3, AlertMsg, TOWMsg_3, TOWWeek3, Top3,
                 Toc3, accuracy3, URAoc_3, URAoc1_3, URAoc2_3, af0_3, af1_3, af2_3 );

   double ClkCorr3 = cc3.svClockBias( dt3 );
   double ClkDrift3 = cc3.svClockDrift( dt3 ); 
   cout << "Clock Bias cc3:  "        << ClkCorr3 << endl;
   cout << "Clock Drift cc3: "        << ClkDrift3 << endl; 
   cout << "Time of Prediction cc3: " << GPSWeekSecond(cc3.getTimeOfPrediction()).printf("%F, %g") << endl;

      // Fourth test case.  Compare against "classic" EngEphemeris
   cout << endl << "Test Case 4: Calculated position using 'classic' EngEphemeris." << endl;
   cout<< "Time= "<< g << endl;
   EngEphemeris EE;
   EE.addSubframe(subframe1, TOWWeek, 3, 1);
   EE.addSubframe(subframe2, TOWWeek, 3, 1);
   EE.addSubframe(subframe3, TOWWeek, 3, 1);

   Xvt xvt = EE.svXvt(dt);
   cout<< "Clock Bias EE:  " << xvt.clkbias << endl;
   cout<< "Clocl Drift EE: " << xvt.clkdrift << endl;

   cout << endl;
   cout << "CC1 Object Dump:" << endl;
   cout << cc1 << endl;

   cout << endl;
   cout << "CC2 Object Dump:" << endl;
   cout << cc2 << endl;

   cout << endl;
   cout << "CC3 Object Dump:" << endl;
   cout << cc3 << endl;

   return(0);
}
