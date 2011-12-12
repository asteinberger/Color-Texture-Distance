/*
  Image.h

  This is a modified version of Gady Agam's code from his
  "Introduction to programming with OpenCV" which was found at:
  http://www.cs.iit.edu/%7Eagam/cs512/lect-notes/opencv-intro/index.html

  Modified by Michael Eckmann 20080219
    - added code to release the IplImage
    - changed b, g and r to blue, green and red
    - changed Rgb to RGB
    - changed Bw to Grey
    - added function to return the imgp getImgPtr()

*/

//#ifndef IMAGE_H_
//#define IMAGE_H_

template<class T> class Image
{
	private:
		IplImage *imgp;

	public:
		Image(IplImage *img=0)
		{
			imgp = img;
		}

		~Image()
		{
			cvReleaseImage(&imgp); // added this to release the image
			imgp = 0;
		}

		void operator = (IplImage *img)
		{
			imgp = img;
		}

		inline T* operator [] (const int rowIndx)
		{
			return ((T*) (imgp->imageData + rowIndx * imgp->widthStep));
		}

		IplImage* getImgPtr()
		{
			return imgp;
		}
};

typedef struct
{
	unsigned char blue, green, red;
} RGBPixel;

typedef struct
{
	float blue, green, red;
} RGBPixelFloat;

typedef Image<RGBPixel>       RGBImage;
typedef Image<RGBPixelFloat>  RGBImageFloat;
typedef Image<unsigned char>  GreyImage;
typedef Image<float>          GreyImageFloat;

//#endif /* IMAGE_H_ */
