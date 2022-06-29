/***********************************************/
/**
* @file digitalFilter.cpp
*
* @brief Digital filter implementation.
*
* @author Matthias Ellmer
* @author Andreas Kvas
* @date 2015-10-29
*
*/
/***********************************************/

#define DOCSTRING_DigitalFilter

#include "base/import.h"
#include "base/fourier.h"
#include "config/configRegister.h"
#include "inputOutput/logging.h"
#include "classes/digitalFilter/digitalFilterMovingAverage.h"
#include "classes/digitalFilter/digitalFilterMedian.h"
#include "classes/digitalFilter/digitalFilterDerivative.h"
#include "classes/digitalFilter/digitalFilterIntegral.h"
#include "classes/digitalFilter/digitalFilterGraceLowpass.h"
#include "classes/digitalFilter/digitalFilterButterworth.h"
#include "classes/digitalFilter/digitalFilterFile.h"
#include "classes/digitalFilter/digitalFilterCorrelation.h"
#include "classes/digitalFilter/digitalFilterWavelet.h"
#include "classes/digitalFilter/digitalFilterNotch.h"
#include "classes/digitalFilter/digitalFilterDecorrelation.h"
#include "classes/digitalFilter/digitalFilterLag.h"
#include "classes/digitalFilter/digitalFilterReduceFilterOutput.h"
#include "classes/digitalFilter/digitalFilter.h"

/***********************************************/

GROOPS_REGISTER_CLASS(DigitalFilter, "digitalFilterType",
                      DigitalFilterMovingAverage,
                      DigitalFilterMedian,
                      DigitalFilterDerivative,
                      DigitalFilterIntegral,
                      DigitalFilterCorrelation,
                      DigitalFilterGraceLowpass,
                      DigitalFilterButterworth,
                      DigitalFilterFile,
                      DigitalFilterWavelet,
                      DigitalFilterNotch,
                      DigitalFilterDecorrelation,
                      DigitalFilterLag,
                      DigitalFilterReduceFilterOutput)

GROOPS_READCONFIG_UNBOUNDED_CLASS(DigitalFilter, "digitalFilterType")

/***********************************************/

DigitalFilter::DigitalFilter(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "create digital filter"))
    {
      if(readConfigChoiceElement(config, "movingAverage", type, "moving average (boxcar) filter"))
        filters.push_back(new DigitalFilterMovingAverage(config));
      if(readConfigChoiceElement(config, "movingMedian",  type, "moving median filter"))
        filters.push_back(new DigitalFilterMedian(config));
      if(readConfigChoiceElement(config, "derivative",    type, "differentiation by polynomial approximation"))
        filters.push_back(new DigitalFilterDerivative(config));
      if(readConfigChoiceElement(config, "integral",      type, "integration by polynomial approximation"))
        filters.push_back(new DigitalFilterIntegral(config));
      if(readConfigChoiceElement(config, "correlation",   type, "correlation by simple coefficient"))
        filters.push_back(new DigitalFilterCorrelation(config));
      if(readConfigChoiceElement(config, "graceLowpass",  type, "GRACE low pass filter (self convolving kernel)"))
        filters.push_back(new DigitalFilterGraceLowpass(config));
      if(readConfigChoiceElement(config, "butterworth",   type, "fixed order digital Butterworth filter"))
        filters.push_back(new DigitalFilterButterworth(config));
      if(readConfigChoiceElement(config, "file",          type, "read ARMA filter from file"))
        filters.push_back(new DigitalFilterFile(config));
      if(readConfigChoiceElement(config, "wavelet",       type, "filter representation of wavelet"))
        filters.push_back(new DigitalFilterWavelet(config));
      if(readConfigChoiceElement(config, "notch",         type, "notch filter"))
        filters.push_back(new DigitalFilterNotch(config));
      if(readConfigChoiceElement(config, "decorrelation",         type, "decorrelation filter"))
        filters.push_back(new DigitalFilterDecorrelation(config));
      if(readConfigChoiceElement(config, "lag",           type, "lag/lead filter"))
        filters.push_back(new DigitalFilterLag(config));
      if(readConfigChoiceElement(config, "reduceFilterOutput",  type, "remove filter output from input signal"))
        filters.push_back(new DigitalFilterReduceFilterOutput(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

DigitalFilter::~DigitalFilter()
{
  for(UInt i=0; i<filters.size(); i++)
    delete filters.at(i);
}

/***********************************************/

Matrix DigitalFilter::filter(const_MatrixSliceRef input) const
{
  try
  {
    Matrix output = input;
    for(UInt i=0; i<filters.size(); i++)
      output = filters.at(i)->filter(output);
    return output;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<std::complex<Double>> DigitalFilter::frequencyResponse(UInt length) const
{
  try
  {
    std::vector<std::complex<Double>> F((length+2)/2, 1.);
    for(UInt i=0; i<filters.size(); i++)
    {
      auto H = filters.at(i)->frequencyResponse(length);
      for(UInt k=0; k<H.size(); k++)
        F.at(k) *= H.at(k);
    }
    return F;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Matrix DigitalFilterARMA::filter(const_MatrixSliceRef input) const
{
  try
  {
    // check input
    // -----------
    if(input.rows()<warmup())
      throw(Exception("Time series is too short (<"+input.rows()%"%i"s+"> elements) to apply a filter with a warmup length of <"+warmup()%"%i"s+">."));

    // filter in frequency domain
    // --------------------------
    if(inFrequencyDomain)
    {
      Matrix padded = pad(input, warmup(), bnStartIndex, padType);
      auto H = frequencyResponse(padded.rows());
      for(UInt k=0; k<padded.columns(); k++) // Filter column-wise
      {
        auto F = Fourier::fft(padded.column(k));
        for(UInt i=0; i<F.size(); i++)
          F.at(i) *= H.at(i);
        copy(Fourier::synthesis(F, (padded.rows()%2==0)), padded.column(k));
      }
      return trim(padded, warmup(), bnStartIndex, padType);
    }

    // filter in time domain
    // ---------------------
    Matrix padded = pad(input, warmup(), bnStartIndex, padType);
    if(backward)
    {
      for(UInt k = 0; k < padded.rows()/2; k++)
        swap(padded.row(k), padded.row(padded.rows() - k - 1));
    }
    Matrix output(padded.rows(), padded.columns());

    const UInt blockSize = std::min(static_cast<UInt>(64), padded.rows());

    Matrix B(bn.rows() + blockSize - 1, blockSize);
    for(UInt k = 0; k < B.columns(); k++)
      copy(bn, B.slice(k, k, bn.rows(), 1));

    for(UInt idxStart = 0; idxStart < output.rows(); idxStart+=blockSize)
    {
      UInt columnCount = std::min(output.rows() - idxStart, blockSize);
      UInt rowCount = std::min(output.rows() - idxStart, B.rows());
      matMult(1.0, B.slice(0, 0, rowCount, columnCount), padded.row(idxStart, columnCount), output.row(idxStart, rowCount));
    }

    if(an.size() > 1)
    {
      Matrix A(an.rows() + blockSize - 1, blockSize);
      for(UInt k = 0; k < A.columns(); k++)
        copy(an, A.slice(k, k, an.rows(), 1));

      MatrixSlice A1(A);
      A1.setType(Matrix::TRIANGULAR, Matrix::LOWER);

      for(UInt idxStart = 0; idxStart < output.rows(); idxStart+=blockSize)
      {
        UInt columnCount = std::min(output.rows() - idxStart, blockSize);
        if(idxStart > 0)
          matMult(-1.0, A.row(blockSize, std::min(an.size() - 1, columnCount)), output.row(idxStart - blockSize, blockSize), output.row(idxStart, std::min(an.size() - 1, columnCount)));

        triangularSolve(1.0, A1.slice(0, 0, columnCount, columnCount), output.row(idxStart, columnCount));
      }
    }

    if(backward)
    {
      for(UInt k = 0; k < output.rows()/2; k++)
        swap(output.row(k), output.row(output.rows() - k - 1));
    }

    return trim(output, warmup(), bnStartIndex, padType);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<std::complex<Double>> DigitalFilterARMA::frequencyResponse(UInt length) const
{
  try
  {
    if((bn.rows() > length) || (an.rows() > length))
      throw(Exception("length must be at least"+std::max(bn.rows(), an.rows())%"%i"s+" for this filter"));

    Vector bPad(length);
    copy(bn.row(bnStartIndex, bn.rows()-bnStartIndex), bPad.row(0, bn.rows()-bnStartIndex));
    copy(bn.row(0, bnStartIndex), bPad.row(bPad.rows()-bnStartIndex, bnStartIndex));

    Vector aPad(length);
    copy(an, aPad.row(0, an.rows()));

    if(backward) // reflect around element one
      for(UInt k=0; k<std::max(bn.rows(), an.rows()); k++)
      {
        std::swap(aPad(k+1), aPad(aPad.rows()-1-k));
        std::swap(bPad(k+1), bPad(bPad.rows()-1-k));
      }

    auto A = Fourier::fft(aPad);
    auto B = Fourier::fft(bPad);
    std::vector<std::complex<Double>> F(A.size());
    for(UInt k=0; k<F.size(); k++)
      F.at(k) = (std::abs(A.at(k))>0.) ? (B.at(k)/A.at(k)) : 1.;

    return F;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt DigitalFilterARMA::warmup() const
{
  return std::max(std::max(bn.rows()-bnStartIndex-1, bnStartIndex), 3*an.rows());
}

/***********************************************/
/***********************************************/

Matrix DigitalFilterBase::pad(const_MatrixSliceRef input, UInt length, UInt timeShift, PadType padType)
{
  try
  {
    if(padType == PadType::NONE)
    {
      if(timeShift > 0)
      {
        Matrix padded(input.rows()+timeShift, input.columns());
        copy(input, padded.row(0, input.rows()));
        return padded;
      }
      else
        return input;
    }

    if(input.rows()<1)
      throw(Exception("Trying to pad a zero length array ("+input.rows()%"%i"s + " x "+ input.columns()%"%i"s+")."));

    Matrix padded(2*length+input.rows()+timeShift, input.columns());
    copy(input, padded.row(length, input.rows()));

    switch(padType)
    {
      case(PadType::ZERO):
        break;

      case(PadType::CONSTANT):
        for(UInt k=0; k<length; k++)
        {
          copy(input.row(0), padded.row(k));
          copy(input.row(input.rows()-1), padded.row(input.rows() + length + k));
        }
        break;

      case(PadType::PERIODIC):
        if(input.rows()<length)
          throw(Exception("Time series is too short ( <"+input.rows()%"%i"s+"> elements) to apply periodic padding for a filter with a warmup length of <"+length%"%i"s+">."));
        for(UInt k=0; k<length; k++)
        {
          copy(input.row(input.rows()-length+k), padded.row(k));
          copy(input.row(k), padded.row(input.rows() + length + k));
        }
        break;

      case(PadType::SYMMETRIC):
        if(input.rows()<length+1)
          throw(Exception("Time series is too short ( <"+input.rows()%"%i"s+"> elements) to apply symmetric padding for a filter with a warmup length of <"+length%"%i"s+">."));
        for(UInt k=0; k<length; k++)
        {
          copy(input.row(k+1), padded.row(length-1-k));
          copy(input.row(input.rows()-2-k), padded.row(input.rows() + length + k));
        }
        break;

      default:
        throw(Exception("Unknown pad type."));
    }

    return padded;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix DigitalFilterBase::trim(const_MatrixSliceRef input, UInt length, UInt timeShift, PadType padType)
{
  try
  {
    if(padType == PadType::NONE)
    {
      if(timeShift > 0)
      {
        return Matrix(input.row(timeShift, input.rows() - timeShift));
      }
      else
        return input;
    }

    return Matrix(input.row(length + timeShift, input.rows()-2*length - timeShift));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, DigitalFilterBase::PadType &padType, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    std::string choice;
    if(readConfigChoice(config, name, choice, mustSet, defaultValue, annotation))
    {
      if(readConfigChoiceElement(config, "none",      choice, "no padding is applied"))                          padType = DigitalFilterBase::PadType::NONE;
      if(readConfigChoiceElement(config, "zero",      choice, "zero padding"))                                   padType = DigitalFilterBase::PadType::ZERO;
      if(readConfigChoiceElement(config, "constant",  choice, "pad using first and last value"))                 padType = DigitalFilterBase::PadType::CONSTANT;
      if(readConfigChoiceElement(config, "periodic",  choice, "periodic continuation of matrix"))                padType = DigitalFilterBase::PadType::PERIODIC;
      if(readConfigChoiceElement(config, "symmetric", choice, "symmetric continuation around the matrix edges")) padType = DigitalFilterBase::PadType::SYMMETRIC;
      endChoice(config);
      return TRUE;
    }

    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
