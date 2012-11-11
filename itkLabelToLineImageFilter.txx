/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLabelToLineImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-16 16:45:10 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLabelToLineImageFilter_txx
#define __itkLabelToLineImageFilter_txx

#include "itkLabelToLineImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkListSample.h"
#include "itkCovarianceSampleFilter.h"
//#include "itkCovarianceCalculator.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkCrossHelper.h"

#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Constructor
 */
template < typename  TInput, typename TOutput  >
LabelToLineImageFilter< TInput, TOutput >
::LabelToLineImageFilter()
{
  this->m_Label = 1;
}


template < typename  TInput, typename TOutput  >
void 
LabelToLineImageFilter< TInput, TOutput >
::GenerateData()
{
  itkDebugMacro(<< "LabelToLineImageFilter generating data ");
  
  typename InputImageType::ConstPointer input = this->GetInput();
  typename OutputImageType::Pointer output = this->GetOutput();

  ImageRegionConstIterator<InputImageType> it;
  it = ImageRegionConstIterator<InputImageType>( input, input->GetRequestedRegion() );
  ImageRegionIterator<OutputImageType> oit;
  this->AllocateOutputs();
  oit = ImageRegionIterator<OutputImageType>(output,
                                             output->GetRequestedRegion());

  output->FillBuffer(static_cast<OutputPixelType>(0));
  
  typedef itk::Statistics::ListSample< VectorType > PointListType;
  typedef std::map< InputPixelType, PointListType::Pointer > PointListMapType;

  //-- Create lists of points for selected label
  PointListType::Pointer plist = PointListType::New();

  oit.GoToBegin();
  it.GoToBegin();
  int found = 0;

  while (!it.IsAtEnd())
    {
    InputPixelType pix = it.Get();
    typename InputImageType::IndexType index = it.GetIndex();

    if (pix == this->m_Label)
      {
      found = 1;
      VectorType mv;
      typename InputImageType::PointType point;
      input->TransformIndexToPhysicalPoint (index, point);
      mv[0] = (double)point[0];
      mv[1] = (double)point[1];
      mv[2] = (double)point[2];
      plist->PushBack(mv);
      oit.Set(it.Get());
      }
    ++it;
    ++oit;
    }

  //-- For each label, perform principal component analysis
  VectorType needleNorm;      // Orientation normal vector of the needle
  VectorType needleTip;      // Tip of the needle closest to the default point.

  InputPixelType pix = this->m_Label;
  PointListType::Pointer sample = plist; 
  std::cout << "=== Label " << pix << "===" << std::endl;
  
  typedef itk::Statistics::CovarianceSampleFilter< PointListType > 
    CovarianceAlgorithmType;
  CovarianceAlgorithmType::Pointer covarianceAlgorithm = 
    CovarianceAlgorithmType::New();
  
  covarianceAlgorithm->SetInput( sample );
  covarianceAlgorithm->Update();
  
  std::cout << "Sample covariance = " << std::endl ; 
  std::cout << covarianceAlgorithm->GetCovarianceMatrix() << std::endl;
  
  CovarianceAlgorithmType::MeasurementVectorType meanVector;
  meanVector = covarianceAlgorithm->GetMean();
  std::cout << "Sample mean = " << meanVector << std::endl ; 
  
  // Perform Symmetric Eigen Analysis
  typedef itk::FixedArray< double, 3 > EigenValuesArrayType;
  typedef itk::Matrix< double, 3, 3 > EigenVectorMatrixType;
  typedef itk::SymmetricEigenAnalysis< CovarianceAlgorithmType::MatrixType,
    EigenValuesArrayType, EigenVectorMatrixType > SymmetricEigenAnalysisType;
  SymmetricEigenAnalysisType analysis ( 3 );
  EigenValuesArrayType eigenValues;
  EigenVectorMatrixType eigenMatrix;
  analysis.SetOrderEigenMagnitudes( true );
  analysis.ComputeEigenValuesAndVectors( covarianceAlgorithm->GetCovarianceMatrix(),
                                         eigenValues, eigenMatrix );    
  
  std::cout << "EigenValues: " << eigenValues << std::endl;
  std::cout << "EigenVectors (each row is an an eigen vector): " << std::endl;
  std::cout << eigenMatrix << std::endl;
  
  // Check the direction of principal component
  VectorType principalVector = eigenMatrix[2];
  double ip = principalVector * m_Normal;
  
  // Calculate the needle orientation vector.
  // If the default vector (m_Normal) and the principal
  // vector are oposit, flip the orientation.
  if (ip >= 0)
    {
    needleNorm = principalVector;
    }
  else
    {
    needleNorm = - principalVector;
    }
  
  needleNorm.Normalize();
  
  PointListType::Iterator iter = sample->Begin();
  
  // To detect the edge of the needle artifact, calculate
  // projections of the points in the needle artifact, and
  // find the farest from the center point (meanVector)
  VectorType vector = iter.GetMeasurementVector();
  double min = (vector-meanVector)*needleNorm;
  double max = min;
  
  while (iter != sample->End())
    {
    vector = iter.GetMeasurementVector();
    double p = (vector-meanVector)*needleNorm;
    if (p < min)
      {
      min = p;
      }
    else if (p > max)
      {
      max = p;
      }
    
    typename InputImageType::PointType point;
    typename OutputImageType::IndexType index;
    point[0] = vector[0];
    point[1] = vector[1];
    point[2] = vector[2];
    output->TransformPhysicalPointToIndex (point, index);
    output->SetPixel(index, pix);
    ++ iter;
    }
  
  needleTip = meanVector + needleNorm * max;
  
  // Output position and orientation of the needle as Affine transform
  // (suppose an identitiy transform when the needle orients (0, 0, 1) direction)
  m_NeedleTransform = NeedleTransformType::New();
  m_NeedleTransform->SetIdentity();

  if (found > 0)
    {
    VectorType nx;
    VectorType ny;
    nx[0] = 1.0; nx[1] = 0.0; nx[2] = 0.0;
    ny[0] = 0.0; ny[1] = 1.0; ny[2] = 0.0;
    typedef itk::CrossHelper< VectorType > CrossType;
    CrossType cross;
    VectorType t = cross(needleNorm, nx);
    VectorType s = cross(t, needleNorm);
    
    t.Normalize();
    s.Normalize();
    needleNorm.Normalize();
    
    NeedleTransformType::MatrixType matrix;
    matrix[0][0] = s[0];
    matrix[0][1] = s[1];
    matrix[0][2] = s[2];
    
    matrix[1][0] = t[0];
    matrix[1][1] = t[1];
    matrix[1][2] = t[2];
    
    matrix[2][0] = needleNorm[0];
    matrix[2][1] = needleNorm[1];
    matrix[2][2] = needleNorm[2];
    
    // Convert translation for slicer coordinate
    m_NeedleTransform->SetMatrix( matrix );
    m_NeedleTransform->Translate( - (matrix * needleTip) );
    }

}


template < typename  TInput, typename TOutput  >
void
LabelToLineImageFilter< TInput, TOutput >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


} // end namespace itk
  
#endif
