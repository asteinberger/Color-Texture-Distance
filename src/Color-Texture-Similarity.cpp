/*
 Morphology.cpp
 written by Adam Steinberger
 */
#include <fstream>
#include <iostream>
#include <sstream>
#include <queue>
#include "cv.h"
#include "highgui.h"
#include "Image.h"

using namespace std;

float L1(float bin1[], float bin2[], int size);
float L2(float bin1[], float bin2[], int size);
int arraySum(int x [], int size);
float arraySum(float x [], int size);
float match(int bin1[], int bin2[]);
int intersection(int bin1[], int bin2[]);
void getLocalBinParts(GreyImage& inputImage, int bins []);
void getColorBins(IplImage *inputImage, int bins []);
GreyImage* crossCorrelate(GreyImage& inputImage, int maskX[3][3], int maskY[3][3],
		float gradBins[], float angleBins[]);

typedef struct {
	int r, c;
} Pixel;

const float pi = 3.14159;
const unsigned char white = 255;
const unsigned char black = 0;
const unsigned char halfWay = 128;

class ScoreCard {
private:
	string filename;
	double score;
public:
	ScoreCard(string f, double s) {
		filename = f;
		score = s;
	} // end ScoreCard constructor
	inline string getFilename() { return filename; }
	inline double getScore() { return score; }
	bool operator() (const ScoreCard& lhs, const ScoreCard&rhs) const {
		return (lhs.score-rhs.score) > 0;
	} // end operator()
	bool operator< (const ScoreCard& card) const {
		return (score<card.score);
	} // end operator()
}; // end ScoreCard class

int main() {

	// setup prewitt mask
	int prewittXmtx [3][3];
	prewittXmtx[0][0] = -1;
	prewittXmtx[0][1] = 0;
	prewittXmtx[0][2] = 1;
	prewittXmtx[1][0] = -1;
	prewittXmtx[1][1] = 0;
	prewittXmtx[1][2] = 1;
	prewittXmtx[2][0] = -1;
	prewittXmtx[2][1] = 0;
	prewittXmtx[2][2] = 1;

	int prewittYmtx [3][3];
	prewittYmtx[0][0] = 1;
	prewittYmtx[0][1] = 1;
	prewittYmtx[0][2] = 1;
	prewittYmtx[1][0] = 0;
	prewittYmtx[1][1] = 0;
	prewittYmtx[1][2] = 0;
	prewittYmtx[2][0] = -1;
	prewittYmtx[2][1] = -1;
	prewittYmtx[2][2] = -1;

	// setup sobel mask
	int sobelXmtx [3][3];
	sobelXmtx[0][0] = -1;
	sobelXmtx[0][1] = 0;
	sobelXmtx[0][2] = 1;
	sobelXmtx[1][0] = -2;
	sobelXmtx[1][1] = 0;
	sobelXmtx[1][2] = 2;
	sobelXmtx[2][0] = -1;
	sobelXmtx[2][1] = 0;
	sobelXmtx[2][2] = 1;

	int sobelYmtx [3][3];
	sobelYmtx[0][0] = 1;
	sobelYmtx[0][1] = 2;
	sobelYmtx[0][2] = 1;
	sobelYmtx[1][0] = 0;
	sobelYmtx[1][1] = 0;
	sobelYmtx[1][2] = 0;
	sobelYmtx[2][0] = -1;
	sobelYmtx[2][1] = -2;
	sobelYmtx[2][2] = -1;

	// open imageNamesFile.txt
	ifstream inFile("imageNamesFile.txt", ios_base::in);
	if (inFile == NULL) {
		cout << "Couldn't open imageNamesFile.txt" << endl;
		return 1;
	} // end if

	// copy all fileNames in txt file to String array
	int size=0;
	string fileNames [100];
	while (!inFile.eof()) {
		inFile >> fileNames[size];
		size++;
	} // end while
	inFile.close();

	// setup bins, load images and run tests
	int *queryBins;
	queryBins = new int[64];
	int *imgBins;
	imgBins = new int[64];
	float *queryGradBins;
	queryGradBins = new float[2];
	float *queryAngleBins;
	queryAngleBins = new float[3];
	int *queryLbpBins;
	queryLbpBins = new int[255];
	float *gradBins;
	gradBins = new float[2];
	float *angleBins;
	angleBins = new float[3];
	int *lbpBins;
	lbpBins = new int[255];
	const char* aWord;
	priority_queue<ScoreCard> scores;

	// load color query image from file
	aWord = fileNames[0].c_str();
	IplImage *queryImg = cvLoadImage(aWord, 3);
	if (!queryImg)
		cout << "Could not load image file: " << aWord << endl;

	// load greyscale query image from file
	IplImage *bwImg = cvLoadImage(aWord, 0);
	if (!bwImg)
			cout << "Could not load image file: " << aWord << endl;
	GreyImage bwImage(bwImg);

	cout << "fileNames[0]: " << aWord << endl;

	// get color histogram bins for image
	getColorBins(queryImg, queryBins);

	// get normalized bins
	float sum = (float) arraySum(queryBins,64);
	float normQueryBins [64];
	for (int i = 0; i < 64; i++) {
		normQueryBins[i] = (float) queryBins[i]/sum;
	} // end for

	// cross-correlate greyscale image with prewitt mask
	GreyImage *crossCor;
	crossCor = crossCorrelate(bwImage,prewittXmtx,prewittYmtx,queryGradBins,queryAngleBins);

	// get normalized bins
	float normQueryGradBins1 [5];
	sum = (float) arraySum(queryGradBins,2) + arraySum(queryAngleBins,3);
	normQueryGradBins1[0] = (float) queryGradBins[0]/sum;
	normQueryGradBins1[1] = (float) queryGradBins[1]/sum;
	normQueryGradBins1[2] = (float) queryAngleBins[0]/sum;
	normQueryGradBins1[3] = (float) queryAngleBins[1]/sum;
	normQueryGradBins1[4] = (float) queryAngleBins[2]/sum;

	// cross-correlate greyscale image with sobel mask
	crossCor = crossCorrelate(bwImage,sobelXmtx,sobelYmtx,queryGradBins,queryAngleBins);

	// get normalized bins
	float normQueryGradBins2 [5];
	sum = (float) arraySum(queryGradBins,2) + arraySum(queryAngleBins,3);
	normQueryGradBins2[0] = (float) queryGradBins[0]/sum;
	normQueryGradBins2[1] = (float) queryGradBins[1]/sum;
	normQueryGradBins2[2] = (float) queryAngleBins[0]/sum;
	normQueryGradBins2[3] = (float) queryAngleBins[1]/sum;
	normQueryGradBins2[4] = (float) queryAngleBins[2]/sum;

	// get local binary partitions for query image
	getLocalBinParts(bwImage,queryLbpBins);

	// get normalized bins
	sum = (float) arraySum(queryLbpBins,255);
	float normQueryLbpBins [255];
	for (int i = 0; i < 255; i++) {
		normQueryLbpBins[i] = (float) queryLbpBins[i]/sum;
	} // end for

	cout << "----------" << endl;

	// compare all other images to query image
	for (int index = 1; index < size; index++) {

		// reset score to 0
		double score = 0;

		// load color image from file
		aWord = fileNames[index].c_str();
		IplImage *inImg = cvLoadImage(aWord, 3);
		if (!inImg)
			cout << "Could not load image file: " << aWord << endl;

		// load greyscale image from file
		IplImage *bwImg = cvLoadImage(aWord, 0);
		if (!bwImg)
				cout << "Could not load image file: " << aWord << endl;

		cout << "fileNames[" << index << "]: " << aWord << endl;

		// get color histogram bins for image
		getColorBins(inImg, imgBins);

		// get normalized bins
		float sum = (float) arraySum(imgBins,64);
		float normBins [64];
		for (int i = 0; i < 64; i++)
			normBins[i] = (float) imgBins[i]/sum;

		// find intersection of image and query image
		int xbin = intersection(queryBins,imgBins);
		cout << "INTERSECTION = " << xbin << endl;

		// add to score
		score += (double) xbin;

		// find match value of image and query image
		float xbin2 = match(queryBins,imgBins);
		cout << "MATCH = " << xbin2 << endl;

		// add to score
		score += (double) xbin2;

		// find l1 distance of color bins from image to query image
		float dist1 = L1(normQueryBins,normBins,64);
		cout << "L1 color = " << dist1 << endl;

		// add to score
		score += (double) dist1;

		// find l2 distance of color bins from image to query image
		float dist2 = L2(normQueryBins,normBins,64);
		cout << "L2 color = " << dist2 << endl;

		// add to score
		score += (double) dist2;

		// cross-correlate greyscale image with prewitt mask
		GreyImage bwImage(bwImg);
		GreyImage *crossCor;
		crossCor = crossCorrelate(bwImage,prewittXmtx,prewittYmtx,gradBins,angleBins);
		string f = fileNames[index] + ".ccp.jpg";
		if (!cvSaveImage(f.c_str(), crossCor->getImgPtr()))
			printf("Couldn't save image.");

		// get normalized bins
		float normGradBins [5];
		sum = (float) arraySum(gradBins,2) + arraySum(angleBins,3);
		normGradBins[0] = (float) gradBins[0]/sum;
		normGradBins[1] = (float) gradBins[1]/sum;
		normGradBins[2] = (float) angleBins[0]/sum;
		normGradBins[3] = (float) angleBins[1]/sum;
		normGradBins[4] = (float) angleBins[2]/sum;

		// find l1 distance of gradient bins from image to query image
		dist1 = L1(normQueryGradBins1,normGradBins,5);
		cout << "L1 prewitt = " << dist1 << endl;

		// add to score
		score += (double) dist1;

		// find l2 distance of gradient bins from image to query image
		dist2 = L2(normQueryGradBins1,normGradBins,5);
		cout << "L2 prewitt = " << dist2 << endl;

		// add to score
		score += (double) dist2;

		// cross-correlate greyscale image with sobel mask
		crossCor = crossCorrelate(bwImage,sobelXmtx,sobelYmtx,gradBins,angleBins);
		string f2 = fileNames[index] + ".ccs.jpg";
		if (!cvSaveImage(f2.c_str(), crossCor->getImgPtr()))
			printf("Couldn't save image.");

		// get normalized bins
		sum = (float) arraySum(gradBins,2) + arraySum(angleBins,3);
		normGradBins[0] = (float) gradBins[0]/sum;
		normGradBins[1] = (float) gradBins[1]/sum;
		normGradBins[2] = (float) angleBins[0]/sum;
		normGradBins[3] = (float) angleBins[1]/sum;
		normGradBins[4] = (float) angleBins[2]/sum;

		// find l1 distance of gradient bins from image to query image
		dist1 = L1(normQueryGradBins2,normGradBins,5);
		cout << "L1 sobel = " << dist1 << endl;

		// add to score
		score += (double) dist1;

		// find l2 distance of gradient bins from image to query image
		dist2 = L2(normQueryGradBins2,normGradBins,5);
		cout << "L2 sobel = " << dist2 << endl;

		// add to score
		score += (double) dist2;

		// get local binary partitions of greyscale image
		getLocalBinParts(bwImage,lbpBins);

		// get normalized bins
		sum = (float) arraySum(lbpBins,255);
		float normLbpBins [255];
		for (int i = 0; i < 255; i++) {
			normLbpBins[i] = (float) lbpBins[i]/sum;
		} // end for

		// find l1 distance of lbp bins from image to query image
		dist1 = L1(normQueryLbpBins,normLbpBins,255);
		cout << "L1 lbp = " << dist1 << endl;

		// add to score
		score += (double) dist1;

		// find l2 distance of lbp bins from image to query image
		dist2 = L2(normQueryLbpBins,normLbpBins,255);
		cout << "L2 lbp = " << dist2 << endl;

		// add to score
		score += (double) dist2;

		cout << "SCORE = " << score << endl;

		// make new scorecard to keep track of image scores
		ScoreCard s = ScoreCard(fileNames[index],score);
		scores.push(s);

		cout << "----------" << endl;

	} // end for

	// print images in order of highest to lowest score
	for (int index = 0; index < size; index++) {
		ScoreCard s = scores.top();
		scores.pop();
		cout << s.getFilename() << " has score " << s.getScore() << endl;
	} // end for

	return (0);

} // end main()

void getLocalBinParts(GreyImage& inputImage, int bins []) {

	// set all elements of bins to zero
	for (int i = 0; i < 255; i++)
		bins[i] = 0;

	// loop through all pixels of image
	for (int row = 0; row < inputImage.getImgPtr()->height; row++) {
		for (int col = 0; col < inputImage.getImgPtr()->width; col++) {

			int parts = 0;

			// get pixel neighbors
			Pixel neighbors [8];
			neighbors[0].r = row-1;
			neighbors[0].c = col-1;
			neighbors[1].r = row-1;
			neighbors[1].c = col;
			neighbors[2].r = row-1;
			neighbors[2].c = col+1;
			neighbors[3].r = row;
			neighbors[3].c = col+1;
			neighbors[4].r = row+1;
			neighbors[4].c = col+1;
			neighbors[5].r = row+1;
			neighbors[5].c = col;
			neighbors[6].r = row+1;
			neighbors[6].c = col-1;
			neighbors[7].r = row;
			neighbors[7].c = col-1;

			// for each neighbor of pixel, binary partition is 1 if neighbor pixel's intensity
			// is greater than current pixel's intensity
			for (int n = 0; n < 8; n++) {
				Pixel p = neighbors[n];
				if (inputImage[p.r][p.c] > inputImage[row][col]) {
					parts += (int) pow(2,n);
				} // end if
			} // end for

			// add one to the correct partition bin
			bins[parts]++;

		} // end for
	} // end for

} // end getLocalBinParts()

float L1(float bin1[], float bin2[], int size) {
	float sum = 0;
	for (int i = 0; i < size; i++) {
		sum += abs(bin1[i]-bin2[i]);
	} // end for
	return sum;
} // end L1()

float L2(float bin1[], float bin2[], int size) {
	float sum = 0;
	for (int i = 0; i < size; i++) {
		sum += pow(bin1[i]-bin2[i],2);
	} // end for
	return sqrt(sum);
} // end L2()

int arraySum(int x [], int size) {
	int sum = 0;
	for (int i = 0; i < size; i++) {
		sum += x[i];
	} // end for
	return sum;
} // end sum()

float arraySum(float x [], int size) {
	float sum = 0;
	for (int i = 0; i < size; i++) {
		sum += x[i];
	} // end for
	return sum;
} // end sum()

int intersection(int bin1[], int bin2[]) {
	int sum = 0;
	for (int i = 0; i < 64; i++) {
		sum += min(bin1[i],bin2[i]);
	} // end for
	return sum;
} // end intersection

float match(int bin1[], int bin2[]) {
	int sum = 0;
	for (int i = 0; i < 64; i++) {
		sum += min(bin1[i],bin2[i]);
	} // end for
	int model = arraySum(bin2,64);
	float normSum = (float) sum / model;
	return normSum;
} // end intersection

void getColorBins(IplImage *inputImage, int bins[]) {

	// set all elements of bins to zero
	for (int i = 0; i < 64; i++)
		bins[i] = 0;

	// loop through all pixels of image
	for (int row = 0; row < inputImage->height; row++) {
		for (int col = 0; col < inputImage->width; col++) {

			// create scalar of current pixel
			CvScalar s;
			s = cvGet2D(inputImage,row,col);

			// get color values of current pixel
			int blue = s.val[0];
			int green = s.val[1];
			int red = s.val[2];

			// get 2 highest bits of blue value
			int highBitBlue = blue/128;
			int temp = blue-highBitBlue*128;
			int lowBitBlue = temp/64;

			// get 2 highest bits of green value
			int highBitGreen = green/128;
			temp = green-highBitGreen*128;
			int lowBitGreen = temp/64;

			// get 2 highest bits of red value
			int highBitRed = red/128;
			temp = red-highBitRed*128;
			int lowBitRed = temp/64;

			// determine correct bin for pixel color
			int bin = 32*highBitBlue + 16*lowBitBlue + 8*highBitGreen + 4*lowBitGreen
					+ 2*highBitRed + lowBitRed;

			// add one to correct color bin
			bins[bin]++;

		} // end for
	} // end for

} // end getColorBins()

GreyImage* crossCorrelate(GreyImage& inputImage, int maskX[3][3], int maskY[3][3], float gradBins[], float angleBins[])  {

	// get input image dimensions
	int inputW = inputImage.getImgPtr()->width;
	int inputH = inputImage.getImgPtr()->height;

	// setup output image dimensions
	CvSize outputImgSize;
	outputImgSize.width = inputW;
	outputImgSize.height = inputH;

	// create output image
	IplImage* outIm = cvCreateImage(outputImgSize, IPL_DEPTH_8U, 1);
	GreyImage* outImage = new GreyImage(outIm);

	// clear out gradient bins
	gradBins[0] = 0;
	gradBins[1] = 0;
	angleBins[0] = 0;
	angleBins[1] = 0;
	angleBins[2] = 0;

	int inImgRow, inImgCol;

	// loop through every pixel of image
	for (int iiRow = 0; iiRow < outputImgSize.height; iiRow++) {

		for (int iiCol = 0; iiCol < outputImgSize.width; iiCol++) {

			// clear out X and Y sums
			float sumX = 0;
			float sumY = 0;

			// loop through 3x3 pixel neighborhood
			for (int row = 0; row < 3; row++) {

				for (int col = 0; col < 3; col++) {

					// get current neighbor coordinates
					inImgRow = iiRow - 1 + row;
					inImgCol = iiCol - 1 + col;

					// as long as neighbor is in image, cross-correlate image with 3x3 mask
					if ((inImgRow >= 0) && (inImgRow < inputH) &&
							(inImgCol >= 0) && (inImgCol < inputW)) {
						sumX += (float) inputImage[inImgRow][inImgCol] * maskX[row][col];
						sumY += (float) inputImage[inImgRow][inImgCol] * maskY[row][col];
					} // end if

				} // end for

			} // end for

			// find change in X and Y
			float deltaX = sumX/6.0;
			float deltaY = sumY/6.0;

			// find gradient magnitude and angle
			float gradient = sqrt(pow(deltaX,2)+pow(deltaY,2));
			float angle = atan2(deltaY,deltaX);

			// set output image pixel to gradient magnitude
			(*outImage)[iiRow][iiCol] = (unsigned char) gradient;

			// add one to correct gradient magnitude bin
			if (gradient < 50) {
				gradBins[0]++;
			} else {
				gradBins[1]++;
			} // end if

			// add one to correct gradient angle bin
			if ((angle <= pi/8.0) && (angle >= -pi/8.0)) {
				angleBins[0]++;
			} else if (((angle < 3.0*pi/8.0) && (angle > pi/8.0))
					|| ((angle > -3.0*pi/8.0) && (angle < -pi/8.0))) {
				angleBins[1]++;
			} else {
				angleBins[2]++;
			} // end if

		} // end for

	} // end for

	return outImage;

} // end crossCorrelate()
