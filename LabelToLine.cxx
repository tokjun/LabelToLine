#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkSymmetricSecondRankTensor.h"
#include "itkLabelToLineImageFilter.h"
#include "itkTransformFileWriter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"


#include "itkPluginUtilities.h"
#include "LabelToLineCLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace {

template<class T> int DoIt( int argc, char * argv[], T )
{
  PARSE_ARGS;

  const     unsigned int        Dimension       = 3;

  typedef   T                   FileInputPixelType;
  typedef   int                 InternalPixelType;
  typedef   int                 OutputPixelType;

  typedef   itk::Image< FileInputPixelType, Dimension > FileInputImageType;
  typedef   itk::Image< InternalPixelType, Dimension >  InternalImageType;
  typedef   itk::Image< OutputPixelType, Dimension >    OutputImageType;

  typedef   itk::ImageFileReader< InternalImageType >  ReaderType;
  typedef   itk::ImageFileWriter< OutputImageType > WriterType;

  typedef   itk::TransformFileWriter TransformWriterType;

  // Line detection filter
  typedef   itk::LabelToLineImageFilter<
  InternalImageType, OutputImageType > NeedleFilterType;

  typename ReaderType::Pointer reader = ReaderType::New();  
  typename WriterType::Pointer writer = WriterType::New();
  typename NeedleFilterType::Pointer needleFilter = NeedleFilterType::New();
  typename TransformWriterType::Pointer needleTransformWriter = TransformWriterType::New();

  reader->SetFileName( inputVolume.c_str() );
  writer->SetFileName( outputVolume.c_str() );

  needleFilter->SetInput( reader->GetOutput() );
  needleFilter->SetLabel( needlelabel );
  needleFilter->SetNormal (static_cast< double >(-normal[0]),
                           static_cast< double >(-normal[1]),
                           static_cast< double >(normal[2]));

  writer->SetInput( needleFilter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject &err)
    {
    std::cerr << err << std::endl;
    return EXIT_FAILURE ;
    }

  typedef typename NeedleFilterType::NeedleTransformType TransformType;
  TransformType::Pointer transform = needleFilter->GetNeedleTransform();
  needleTransformWriter->SetFileName( needleTransform );
  needleTransformWriter->SetInput( transform );

  try
    {
    needleTransformWriter->Update();
    }
  catch (itk::ExceptionObject &err)
    {
    std::cerr << err << std::endl;
    return EXIT_FAILURE ;
    }

  return EXIT_SUCCESS;
}

} // end of anonymous namespace


int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType (inputVolume, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch (componentType)
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0));
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0));
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0));
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0));
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0));
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0));
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0));
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0));
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0));
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0));
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }

  catch( itk::ExceptionObject &excep)
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
