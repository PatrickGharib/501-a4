#include <ctime>

#include <array>

#include <iostream>

#include <string>

#include <fstream>
 //#include <vector> 
#include <cstdint>

#include <algorithm>

#include <string.h>

#include <math.h>

#include <cmath>

#include <limits.h>

#include <stdio.h>

#include <stdlib.h>

using namespace std;

//struct to organize wave info from input 
typedef struct WAVEHEADER {
  //TODO change to int and stuff for these cause its apparently really bad 

  char RIFF[4];
  int ChunkSize;
  char format[4];
  char Subchunk1ID[4];
  int Subchunk1Size;
  short AudioFormat;
  short NumChannels;
  int SampleRate;
  int ByteRate;
  short BlockAlign;
  short BitsPerSample;

}
wavHdr;
typedef struct DATACHUNK {
  char Subchunk2ID[4];
  int Subchunk2Size;
}
dChunk;
#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr 

//Prototyping
void four1(double data[], int nn, int isign);
void convolve(double* wavDatFreq1, double *wavDatFreq2, double* outFreq, int arrayLength);
size_t fwriteILSB(int data, FILE * fs);
size_t fwriteSLSB(short data, FILE * fs);
int getWaveFileSize(FILE * inFile);
int8_t * getWavData(FILE * wavFile);
void outPutScaler(double * outputArray, short *sig, int outputSize);

int main(int argc, char ** argv) {

  //check user input for right input
  if (argc != 4) {
    cout << "Usage: convolve <inputfile(.wav)> <IRfile(.wav)> <outputFile(.wav)>";
    exit(0);
  }

  wavHdr wavHeader1;
  wavHdr wavHeader2;

  int sizeOfHeader1 = sizeof(wavHeader1);
  int sizeOfHeader2 = sizeof(wavHeader2);

  const char * inputFile = argv[1];
//-------------------------------------------------------------------------------------------------------
 
  FILE * wavFile = fopen(inputFile, "rb");
  if (wavFile == NULL) {
    fprintf(stderr, "Unable to open inputFile: %s\n", inputFile);
    return 1;
  }

  //read in the header
  fread( & wavHeader1, sizeOfHeader1, 1, wavFile);
  if (wavHeader1.Subchunk1Size == 18) {
    short emptyBytes;
    fread((char * ) & emptyBytes, sizeof(short), 1, wavFile);
  }

  dChunk dataChunk1;
  while (true) {
    fread( & dataChunk1, sizeof(dataChunk1), 1, wavFile);

    if ( * (unsigned int * ) & dataChunk1.Subchunk2ID == 0x61746164)
      break;
    //skip chunk data bytes
    fseek(wavFile, dataChunk1.Subchunk2Size, SEEK_CUR);
  }

  char * value = new char[dataChunk1.Subchunk2Size];
  
  fread(value, dataChunk1.Subchunk2Size, 1, wavFile);

  fclose(wavFile);
  short * sig1;
  int sigSize1;
  sig1 =NULL;
  if (wavHeader1.BitsPerSample == 8) {
    sigSize1 = dataChunk1.Subchunk2Size;
    sig1 = new short[sigSize1];
    for (int i = 0; i < dataChunk1.Subchunk2Size; i++) {
      sig1[i] = (short)((unsigned char) value[i]);
    }
  } else {
    sigSize1 = dataChunk1.Subchunk2Size / 2;
    sig1 = new short[sigSize1];
    
    short shortData;
    for (int i = 0; i < dataChunk1.Subchunk2Size; i += 2) {
      shortData = (short)((unsigned char) value[i]);
      shortData = (short)((unsigned char) value[i + 1]) * 256; 
      sig1[i / 2] = shortData;
    }
  }
  
  double * wavData1 = new double[dataChunk1.Subchunk2Size];

  for (int i = 0; i < sigSize1; i++) {
    wavData1[i] = ((double) sig1[i]) / 32678.0;
   }
 //-------------------------------------------------------------------------------------

  const char * iRFile = argv[2];
  wavFile = fopen(iRFile, "rb");
  if (wavFile == NULL) {
    fprintf(stderr, "Unable to open inputFile: %s\n", iRFile);
    return 1;
  }

  //read in the header
  fread( & wavHeader2, sizeOfHeader2, 1, wavFile);
  if (wavHeader2.Subchunk1Size == 18) {
    short emptyBytes;
    fread((char * ) & emptyBytes, sizeof(short), 1, wavFile);
  }

  dChunk dataChunk2;
  while (true) {
    fread( & dataChunk2, sizeof(dataChunk2), 1, wavFile);
    if ( * (unsigned int * ) & dataChunk2.Subchunk2ID == 0x61746164)
      break;
    //skip chunk data bytes
    fseek(wavFile, dataChunk2.Subchunk2Size, SEEK_CUR);
  }

  
  char * value2 = new char[dataChunk2.Subchunk2Size];
  fread(value2, dataChunk2.Subchunk2Size, 1, wavFile);
  fclose(wavFile);

  short * sig2;
  int sigSize2;
  sig2=NULL;
  if (wavHeader2.BitsPerSample == 8) {
    sigSize2 = dataChunk2.Subchunk2Size;
    sig2 = new short[sigSize2];
    for (int i = 0; i < dataChunk2.Subchunk2Size; i++) {
      sig2[i] = (short)((unsigned char) value2[i]);
    }
  } else {
    sigSize2 = dataChunk2.Subchunk2Size / 2;
    sig2 = new short[sigSize2];
    
    short shortData2;
    for (int i = 0; i < dataChunk2.Subchunk2Size; i += 2) {
      shortData2 = (short)((unsigned char) value2[i]);
      shortData2 = (short)((unsigned char) value2[i + 1]) * 256; 
      sig2[i / 2] = shortData2;
    }
  }
  double * wavData2 = new double[dataChunk2.Subchunk2Size];
  for (int i = 0; i < sigSize2; i++) {
    wavData2[i] = ((double) sig2[i]) / 32678.0;
  } 
//---------------------------------------------------------------------------------------------------------------
int mxLength = 0; 
if (sigSize1<=sigSize2) {mxLength = sigSize2;}
else {mxLength = sigSize1;}

int largestPowTwo = 1;
while(largestPowTwo < mxLength){
  largestPowTwo = largestPowTwo << 1;
}

int mxLength2 = largestPowTwo << 1;
double* wavDatFreq1 = new double[mxLength2];
double* wavDatFreq2 = new double[mxLength2];


memset(wavDatFreq1,0.0,mxLength);
memset(wavDatFreq2,0.0,mxLength2);


for(int i = 0; i < sigSize1; i++){wavDatFreq1[i*2] = wavData1[i];}
for(int i = 0; i < sigSize2; i++){wavDatFreq2[i*2] = wavData2[i];}

four1(wavDatFreq1-1,largestPowTwo,1);
four1(wavDatFreq2-1,largestPowTwo,1);


//---------------------------------------------------------------------------------------------------------------
  int outputNumberOfSamples = sigSize1 + sigSize2 - 1;
  
  double * outputArrayFreq = new double[mxLength2];

  convolve(wavDatFreq1, wavDatFreq2, outputArrayFreq, mxLength2);
  four1(outputArrayFreq-1,largestPowTwo,-1);

  double * outputfreq = new double[outputNumberOfSamples];
  double outmx = 0.0;

            for(int i = 0; i < outputNumberOfSamples; i++){
                double sigSamp = outputArrayFreq[i*2];
                outputfreq[i] = sigSamp;
 
                if(outmx < abs(sigSamp)){
                    outmx = abs(sigSamp);
                }
            }
  outPutScaler(outputfreq, sig1, outputNumberOfSamples);
 
//----------------------------------------------------------------------------------------------------
  wavHdr outputWaveHeader;
  dChunk outputChunk;

  outputWaveHeader.ChunkSize = 36 + outputChunk.Subchunk2Size;
  outputWaveHeader.BlockAlign = wavHeader1.NumChannels * (wavHeader1.BitsPerSample / 8);
  outputWaveHeader.ByteRate = (int) wavHeader1.SampleRate * outputWaveHeader.BlockAlign;
  outputChunk.Subchunk2Size = wavHeader1.NumChannels * outputNumberOfSamples * (wavHeader1.BitsPerSample / 8);
  

  const char * outputFile = argv[3];
  wavFile = fopen(outputFile, "wb");


  fputs("RIFF", wavFile);
  fwriteILSB(outputWaveHeader.ChunkSize, wavFile);
  fputs("WAVE", wavFile);
  fputs("fmt ", wavFile);
  fwriteILSB(16, wavFile); 
  fwriteSLSB(1, wavFile); 
  fwriteSLSB(wavHeader1.NumChannels, wavFile);
  fwriteILSB(wavHeader1.SampleRate, wavFile);
  fwriteILSB(outputWaveHeader.ByteRate, wavFile);
  fwriteSLSB(outputWaveHeader.BlockAlign, wavFile);
  fwriteSLSB(wavHeader1.BitsPerSample, wavFile);
  fputs("data", wavFile);
  fwriteILSB(outputChunk.Subchunk2Size, wavFile);
  for (int i = 0; i < outputChunk.Subchunk2Size; i++) {
    fwriteSLSB(outputfreq[i], wavFile);
  }
  
  
  fclose(wavFile);
  return 0;

}
//-------------------------------------------------------------------------------------------------------------------------------------
int getWaveFileSize(FILE * inFile) {
  int fileSize = 0;
  fseek(inFile, 0, SEEK_END);
  fileSize = ftell(inFile);
  fseek(inFile, 0, SEEK_SET);
  return fileSize;
}
//------------------------------------------------------------------------------------------------------------------
void convolve(double* wavDatFreq1, double *wavDatFreq2, double* outFreq, int arrayLength){
    printf("Starting convolution loops...\n");
    
    for(int i = 0; i < arrayLength; i+= 2){
       double real1 = wavDatFreq1[i]; 
       double real2 = wavDatFreq2[i];
       double imaginary1 = wavDatFreq1[i+1];
       double imaginary2 = wavDatFreq2[i+1];
        outFreq[i] = real1 * real2 - imaginary1 * imaginary2 ; 
        outFreq[i+1] = imaginary1 * real2 + real1 * imaginary2; 
     
    }
}
//-----------------------------------------------------------------------------------------------

void four1(double data[], int nn, int isign){
 
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
 
    n = nn << 1;
    j = 1;
 
    for (i = 1; i < n; i += 2) {
    if (j > i) {
        SWAP(data[j], data[i]);
        SWAP(data[j+1], data[i+1]);
    }
    m = nn;
    while (m >= 2 && j > m) {
        j -= m;
        m >>= 1;
    }
    j += m;
    }
 
    mmax = 2;
    while (n > mmax) {
    istep = mmax << 1;
    theta = isign * (6.28318530717959 / mmax);
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1; m < mmax; m += 2) {
        for (i = m; i <= n; i += istep) {
        j = i + mmax;
        tempr = wr * data[j] - wi * data[j+1];
        tempi = wr * data[j+1] + wi * data[j];
        data[j] = data[i] - tempr;
        data[j+1] = data[i+1] - tempi;
        data[i] += tempr;
        data[i+1] += tempi;
        }
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
    }
}
//------------------------------------------------------------------------------------------------------------------
size_t fwriteILSB(int data, FILE * fs) {

  unsigned char charArray[4];
  charArray[3] = (unsigned char)((data >> 24) & 0xFF);
  charArray[2] = (unsigned char)((data >> 16) & 0xFF);
  charArray[1] = (unsigned char)((data >> 8) & 0xFF);
  charArray[0] = (unsigned char)(data & 0xFF);
  return fwrite(charArray, sizeof(unsigned char), 4, fs);
}

//---------------------------------------------------------------------------------------------------------------------------------------
size_t fwriteSLSB(short data, FILE * fs) {
  unsigned char charArray[2];
  charArray[1] = (unsigned char)((data >> 8) & 0xFF);
  charArray[0] = (unsigned char)(data & 0xFF);
  return fwrite(charArray, sizeof(unsigned char), 2, fs);
}
//------------------------------------------------------------------------------------------------------------------------------------- 
void outPutScaler(double* outputArray, short* sig, int outputSize) {
  double maxIn = 0.0;
  double maxOut = 0.0;
  for (int i = 0; i < outputSize; i++) {
    if (sig[i] > maxIn)
      maxIn = sig[i];
    if (outputArray[i] > maxOut)
      maxOut = outputArray[i];
  }
  double inOutRatio = maxin/maxout;
  for (int i = 0; i < outputSize; i++) {
    outputArray[i] = outputArray[i]*inOutRatio;
  }
}