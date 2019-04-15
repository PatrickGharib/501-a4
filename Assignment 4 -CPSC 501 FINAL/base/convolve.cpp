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

//Prototyping
void convolve(double x[], int N, double h[], int M, double y[], int P);
size_t fwriteILSB(int data, FILE * fileStream);
size_t fwriteSLSB(short data, FILE * fileStream);
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
    cout<<sigSize1<<endl;
    sig1 = new short[sigSize1];
    for (int i = 0; i < dataChunk1.Subchunk2Size; i++) {
      sig1[i] = (short)((unsigned char) value[i]);
    }
  } else {
    sigSize1 = dataChunk1.Subchunk2Size / 2;
    cout<<sigSize1<<endl;
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
  int outputNumberOfSamples = sigSize1 + sigSize2 - 1;
  short * outputsig = new short[outputNumberOfSamples];
  double * outputArray = new double[outputNumberOfSamples];

  convolve(wavData1, sigSize1, wavData2, sigSize2, outputArray, outputNumberOfSamples);
    
  outPutScaler(outputArray, sig1, outputNumberOfSamples);

  for (int i = 0; i < outputNumberOfSamples; i++) {
    outputsig[i] = (short) outputArray[i];
  }

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
    fwriteSLSB(outputsig[i], wavFile);
  }
  
  
  fclose(wavFile);
  return 0;

}

int getWaveFileSize(FILE * inFile) {
  int fileSize = 0;
  fseek(inFile, 0, SEEK_END);
  fileSize = ftell(inFile);
  fseek(inFile, 0, SEEK_SET);
  return fileSize;
}

/*****************************************************************************
 *
 *    Function:     convolve
 *    author: Leonard Manzara
 *
 *    Description:  Convolves two signals, producing an output signal.
 *                  The convolution is done in the time domain using the
 *                  "Input Side Algorithm" (see Smith, p. 112-115).
 *
 *    Parameters:   x[] is the signal to be convolved
 *                  N is the number of samples in the vector x[]
 *                  h[] is the impulse response, which is convolved with x[]
 *                  M is the number of samples in the vector h[]
 *                  y[] is the output signal, the result of the convolution
 *                  P is the number of samples in the vector y[].  P must
 *                       equal N + M - 1
 *
 *****************************************************************************/

void convolve(double x[], int N, double h[], int M, double y[], int P) {
  int n, m;

  /*  Make sure the output buffer is the right size: P = N + M - 1  */
  if (P != (N + M - 1)) {
    printf("Output signal vector is the wrong size\n");
    printf("It is %-d, but should be %-d\n", P, (N + M - 1));
    printf("Aborting convolution\n");

    return;
  }
  /*  Clear the output buffer y[] to all zero values  */
  for (n = 0; n < P; n++) {
    y[n] = 0.0;
    //cout << y[n];
  }
  int start;
  int stop;
  bool flag = true;
  /*  Do the convolution  */
  /*  Outer loop:  process each input value x[n] in turn  */
  for (n = 0; n < N; n++) {
    if (flag) {
      start = clock();
      flag = false;
    }
    /*  Inner loop:  process x[n] with each sample of h[]  */
    for (m = 0; m < M; m++) {
      y[n + m] += x[n] * h[m];

      //cout << y[n+m]<<endl;
    //  cout << x[n] <<endl;
    }

    //calculate estimated time to completion
    if (n % 10000 == 0) {
      stop = clock();
      cout << "estimated time : " << (((N - n) / 10000) * ((stop - start) / double(CLOCKS_PER_SEC))) / 60 << endl;
      flag = true;
    }
  }
}


size_t fwriteILSB(int data, FILE * fileStream) {

  unsigned char charArray[4];
  charArray[3] = (unsigned char)((data >> 24) & 0xFF);
  charArray[2] = (unsigned char)((data >> 16) & 0xFF);
  charArray[1] = (unsigned char)((data >> 8) & 0xFF);
  charArray[0] = (unsigned char)(data & 0xFF);
  return fwrite(charArray, sizeof(unsigned char), 4, fileStream);
}


size_t fwriteSLSB(short data, FILE * fileStream) {
  unsigned char charArray[2];
  charArray[1] = (unsigned char)((data >> 8) & 0xFF);
  charArray[0] = (unsigned char)(data & 0xFF);
  return fwrite(charArray, sizeof(unsigned char), 2, fileStream);
}

void outPutScaler(double* outputArray, short* sig, int outputSize) {
  double maxIn = 0.0;
  double maxOut = 0.0;
  for (int i = 0; i < outputSize; i++) {
    if (sig[i] > maxIn)
      maxIn = sig[i];
    if (outputArray[i] > maxOut)
      maxOut = outputArray[i];
  }
  for (int i = 0; i < outputSize; i++) {
    outputArray[i] = outputArray[i] / maxOut * maxIn;
  }
}