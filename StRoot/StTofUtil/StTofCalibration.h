/*******************************************************************
 *
 * $Id: StTofCalibration.h,v 1.1 2001/09/28 19:09:40 llope Exp $
 *
 * Author: Frank Geurts
 *****************************************************************
 *
 * Description: Calibration class for TOFp
 *
 *****************************************************************
 *
 * $Log: StTofCalibration.h,v $
 * Revision 1.1  2001/09/28 19:09:40  llope
 * first version
 *
 *******************************************************************/
#ifndef STTOFCALIBRATION_H
#define STTOFCALIBRATION_H
#include <vector>
#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

#define TOFP_MAX_SLATS 6000

struct StructSlatCalib {
  //   StructSlatCalib() { }
  float offset_tdc;
  float cc_tdc;
  float ods_tdc;
  float offset_adc;
  float cc_adc;
  float ods_adc;
};

#ifndef ST_NO_TEMPLATE_DEF_ARGS
typedef vector<StructSlatCalib> slatCalibVector;
#else
typedef vector<StructSlatCalib, allocator<StructSlatCalib> > slatCalibVector;
#endif   
typedef slatCalibVector::iterator slatCalibIter;


class StTofCalibration{
 private:
  slatCalibVector mSlatCalibVec; //!
 public:
  StTofCalibration();
  ~StTofCalibration();
  void init();
  void print(); 
  int numberOfEntries() const;
  StructSlatCalib slat(int) const;
};

inline StructSlatCalib StTofCalibration::slat(int index) const {return mSlatCalibVec[index];}
inline int StTofCalibration::numberOfEntries() const {return mSlatCalibVec.size();}
#endif
