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
size_t fwriteIntLSB(int data, FILE * fileStream);
size_t fwriteShortLSB(short data, FILE * fileStream);
int getWaveFileSize(FILE * inFile);
int8_t * getWavData(FILE * wavFile);
void scaleOutputSignal(double * outputArray, short *signal, int outputSize);

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
  //TODO need to move this to save some mem so that is not just doubleing
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

    // cout << dataChunk1.Subchunk2ID[0]<< dataChunk1.Subchunk2ID[1]<< dataChunk1.Subchunk2ID[2]<< dataChunk1.Subchunk2ID[3]<<endl;
    //printf("SubChunk2ID: %c\n", dataChunk1.Subchunk2ID[0]);
    printf("%c%c%c%c\t"
      "%li\n", dataChunk1.Subchunk2ID[0], dataChunk1.Subchunk2ID[1], dataChunk1.Subchunk2ID[2], dataChunk1.Subchunk2ID[3], dataChunk1.Subchunk2Size);
    if ( * (unsigned int * ) & dataChunk1.Subchunk2ID == 0x61746164)
      break;
    //skip chunk data bytes
    fseek(wavFile, dataChunk1.Subchunk2Size, SEEK_CUR);
  }

  //Number of samples
  int sample_size = wavHeader1.BitsPerSample / 8;
  int samples_count = dataChunk1.Subchunk2Size * 8 / wavHeader1.BitsPerSample;
  printf("Samples count = %i\n", samples_count);

  char * value = new char[dataChunk1.Subchunk2Size];
  //memset(value, 0, sizeof(char) * samples_count);
  fread(value, dataChunk1.Subchunk2Size, 1, wavFile);
  // for(int i = 0; i < dataChunk1.Subchunk2Size; i++ ){
  //     cout << value[i];
  // }
  fclose(wavFile);
  short * signal1;
  int signalSize1;
  signal1 =NULL;
  if (wavHeader1.BitsPerSample == 8) {
    signalSize1 = dataChunk1.Subchunk2Size;
    cout<<signalSize1<<endl;
    signal1 = new short[signalSize1];
    for (int i = 0; i < dataChunk1.Subchunk2Size; i++) {
      signal1[i] = (short)((unsigned char) value[i]);
    }
  } else {
    signalSize1 = dataChunk1.Subchunk2Size / 2;
    cout<<signalSize1<<endl;
    signal1 = new short[signalSize1];
    //every 2 chars (i..e. a short) is one s siple
    short shortData;
    for (int i = 0; i < dataChunk1.Subchunk2Size; i += 2) {
      shortData = (short)((unsigned char) value[i]);
      shortData = (short)((unsigned char) value[i + 1]) * 256; //shift to next 8 bits
      signal1[i / 2] = shortData;
    }
  }
  //TODO make this into a method as a refactoring
  double * wavData1 = new double[dataChunk1.Subchunk2Size];

  for (int i = 0; i < signalSize1; i++) {
    wavData1[i] = ((double) signal1[i]) / 32678.0;
   //cout<<signal1[i]  ;


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

    // cout << dataChunk1.Subchunk2ID[0]<< dataChunk1.Subchunk2ID[1]<< dataChunk1.Subchunk2ID[2]<< dataChunk1.Subchunk2ID[3]<<endl;
    //printf("SubChunk2ID: %c\n", dataChunk1.Subchunk2ID[0]);
    printf("%c%c%c%c\t"
      "%li\n", dataChunk2.Subchunk2ID[0], dataChunk2.Subchunk2ID[1], dataChunk2.Subchunk2ID[2], dataChunk2.Subchunk2ID[3], dataChunk2.Subchunk2Size);
    if ( * (unsigned int * ) & dataChunk2.Subchunk2ID == 0x61746164)
      break;
    //skip chunk data bytes
    fseek(wavFile, dataChunk2.Subchunk2Size, SEEK_CUR);
  }

  //Number of samples
  int sample_size2 = wavHeader2.BitsPerSample / 8;
  int samples_count2 = dataChunk2.Subchunk2Size * 8 / wavHeader2.BitsPerSample;
  printf("Samples count = %i\n", samples_count2);

  char * value2 = new char[dataChunk2.Subchunk2Size];
  //memset(value, 0, sizeof(char) * samples_count);
  fread(value2, dataChunk2.Subchunk2Size, 1, wavFile);
  // for(int i = 0; i < dataChunk2.Subchunk2Size; i++ ){
  //     cout << value2[i];
  // }
  fclose(wavFile);

  short * signal2;
  int signalSize2;
  signal2=NULL;
  if (wavHeader2.BitsPerSample == 8) {
    signalSize2 = dataChunk2.Subchunk2Size;
    signal2 = new short[signalSize2];
    for (int i = 0; i < dataChunk2.Subchunk2Size; i++) {
      signal2[i] = (short)((unsigned char) value2[i]);
    }
  } else {
    signalSize2 = dataChunk2.Subchunk2Size / 2;
    signal2 = new short[signalSize2];
    //every 2 chars (i.e. a short) is one signal sample
    short shortData2;
    for (int i = 0; i < dataChunk2.Subchunk2Size; i += 2) {
      shortData2 = (short)((unsigned char) value2[i]);
      shortData2 = (short)((unsigned char) value2[i + 1]) * 256; //shift to next 8 bits
      signal2[i / 2] = shortData2;
    }
  }
  double * wavData2 = new double[dataChunk2.Subchunk2Size];

  
//---------------------------------------------------------------------------------------------------------------
  int outputNumberOfSamples = signalSize1 + signalSize2 - 1;
  short * outputSignal = new short[outputNumberOfSamples];
  double * outputArray = new double[outputNumberOfSamples];



  convolve(wavData1, signalSize1, wavData2, signalSize2, outputArray, outputNumberOfSamples);
    
  scaleOutputSignal(outputArray, signal1, outputNumberOfSamples);

  for (int i = 0; i < outputNumberOfSamples; i++) {
    outputSignal[i] = (short) outputArray[i];
  }

  wavHdr outputWaveHeader;
  dChunk outputChunk;
  outputWaveHeader.ChunkSize = 36 + outputChunk.Subchunk2Size;
  outputChunk.Subchunk2Size = wavHeader1.NumChannels * outputNumberOfSamples * (wavHeader1.BitsPerSample / 8);
 

  outputWaveHeader.NumChannels = wavHeader1.NumChannels;
  
  outputWaveHeader.BlockAlign = wavHeader1.NumChannels * (wavHeader1.BitsPerSample / 8);
  outputWaveHeader.ByteRate =(int) wavHeader1.SampleRate * outputWaveHeader.BlockAlign;
  outputWaveHeader.BitsPerSample = wavHeader1.BitsPerSample;
  outputWaveHeader.SampleRate = wavHeader1.SampleRate;
  //int sizeOfHeader3 = sizeof(outputWaveHeader);
  int OUTSIZE = outputChunk.Subchunk2Size;

  const char * outputFile = argv[3];
  wavFile = fopen(outputFile, "wb");
  fputs("RIFF", wavFile);
  fwriteIntLSB(outputWaveHeader.ChunkSize, wavFile);
  fputs("WAVE", wavFile);

  //fmt subchunknew vs malloc
  fputs("fmt ", wavFile);
  fwriteIntLSB(16, wavFile); //subchunk1size should be fixed 16 bytes
  fwriteShortLSB(1, wavFile); // AudioFormat = 1 for PCM
  fwriteShortLSB(outputWaveHeader.NumChannels, wavFile);
  fwriteIntLSB(outputWaveHeader.SampleRate, wavFile);
  fwriteIntLSB(outputWaveHeader.ByteRate, wavFile);
  fwriteShortLSB(outputWaveHeader.BlockAlign, wavFile);
  fwriteShortLSB(outputWaveHeader.BitsPerSample, wavFile);

  //data subchunk
  fputs("data", wavFile);
  fwriteIntLSB(outputChunk.Subchunk2Size, wavFile);
  for (int i = 0; i < OUTSIZE; i++) {
    fwriteShortLSB(outputSignal[i], wavFile);
  }
  
  
  printf("File Type: %s\n", outputWaveHeader.RIFF);
  printf("File Size: %ld\n", outputWaveHeader.ChunkSize);
  printf("WAV Marker: %s\n", outputWaveHeader.format);
  printf("Format Name: %s\n", outputWaveHeader.Subchunk1ID);
  printf("Format Length: %ld\n", outputWaveHeader.Subchunk1Size );
  printf("Format Type: %hd\n", outputWaveHeader.AudioFormat);
  printf("Number of Channels: %hd\n", outputWaveHeader.NumChannels);
  printf("Sample Rate: %ld\n", outputWaveHeader.SampleRate);
  printf("Sample Rate * Bits/Sample * Channels / 8: %ld\n", outputWaveHeader.ByteRate);
  printf("Bits per Sample * Channels / 8.1: %hd\n", outputWaveHeader.BlockAlign);
  printf("Bits per Sample: %hd\n", outputWaveHeader.BitsPerSample);// fwrite(out,outputWaveHeader.BitsPerSample/8,OUTSIZE,wavFile);
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


size_t fwriteIntLSB(int data, FILE * fileStream) {

  unsigned char charArray[4];

  //write int (4 bytes) into fileStream in little-endian
  //little endian writes from least significant byte (LSB) first
  charArray[3] = (unsigned char)((data >> 24) & 0xFF);
  charArray[2] = (unsigned char)((data >> 16) & 0xFF);
  charArray[1] = (unsigned char)((data >> 8) & 0xFF);
  charArray[0] = (unsigned char)(data & 0xFF);

  //use charArray to write values as characters to file
  return fwrite(charArray, sizeof(unsigned char), 4, fileStream);
}


size_t fwriteShortLSB(short data, FILE * fileStream) {

  unsigned char charArray[2];

  //write short (2 bytes) into fileStream in little-endian
  //little endian writes from least significant byte (LSB) first
  charArray[1] = (unsigned char)((data >> 8) & 0xFF);
  charArray[0] = (unsigned char)(data & 0xFF);

  //use charArray to write values as characters to file
  return fwrite(charArray, sizeof(unsigned char), 2, fileStream);
}

void scaleOutputSignal(double* outputArray, short* signal, int outputSize) {
  double inputMaxValue = 0.0;
  double outputMaxValue = 0.0;

  //check for max value in both original and output signals
  for (int i = 0; i < outputSize; i++) {
    
    if (signal[i] > inputMaxValue)
      inputMaxValue = signal[i];

    if (outputArray[i] > outputMaxValue)
      outputMaxValue = outputArray[i];
  }

  for (int i = 0; i < outputSize; i++) {
    outputArray[i] = outputArray[i] / outputMaxValue * inputMaxValue;
  }
}
void shortToDoubleConversion(int signalSize, short signal[]{
  for (int i = 0; i < signalSize; i++) {
    wavData2[i] = ((double) signal[i]) / 32678.0;
  }
}
/*
 void readWavData(wavHdr &wavHeader1, int sizeOfHeader, double* wavData1, const char* fileName){
   FILE* wavFile = fopen(fileName, "r");

    if (wavFile == NULL){
        fprintf(stderr, "Unable to open inputFile: %s\n", fileName);
        return;
    }
  fread(&wavHeader1, 1, sizeOfHeader, wavFile);
  static const uint64_t BUFFERSIZE = wavHeader1.Subchunk2Size;
    char* buffer = new char[BUFFERSIZE];
    wavData1 = new double[BUFFERSIZE];
    fread(buffer, sizeof buffer[0], BUFFERSIZE / (sizeof buffer[0]), wavFile);
        //TODO make this into a method as a refactoring
    for(int i = 0; i < BUFFERSIZE; i++){wavData1[i] = (double) buffer[i]/32678.0;}
    delete[] buffer; 
    buffer = nullptr;

    fclose(wavFile);
}*/
 /*
  printf("WAV File Header read:\n");
  printf("File Type: %s\n", wavHeader1.RIFF);
  printf("File Size: %ld\n", wavHeader1.ChunkSize);
  printf("WAV Marker: %s\n", wavHeader1.format);)
  printf("Format Name: %s\n", wavHeader1.Subchunk1ID);
  printf("Format Length: %ld\n", wavHeader1.Subchunk1Size );
  printf("Format Type: %hd\n", wavHeader1.AudioFormat);
  printf("Number of Channels: %hd\n", wavHeader1.NumChannels);
  printf("Sample Rate: %ld\n", wavHeader1.SampleRate);
  printf("Sample Rate * Bits/Sample * Channels / 8: %ld\n", wavHeader1.ByteRate);
  printf("Bits per Sample * Channels / 8.1: %hd\n", wavHeader1.BlockAlign);
  printf("Bits per Sample: %hd\n", wavHeader1.BitsPerSample);
  */

  ///delete [] buffer1; 
   // outputWaveHeader.RIFF[0] =  'R';
  // outputWaveHeader.RIFF[1] =  'I';
  // outputWaveHeader.RIFF[2] =  'F';
  // outputWaveHeader.RIFF[3] =  'F';

  // outputWaveHeader.format[0] = 'W';
  // outputWaveHeader.format[1] = 'A';
  // outputWaveHeader.format[2] = 'V';
  // outputWaveHeader.format[3] = 'E';

  // outputWaveHeader.Subchunk1ID[0] = 'f';
  // outputWaveHeader.Subchunk1ID[1] = 'm';
  // outputWaveHeader.Subchunk1ID[2] = 't';

  // outputChunk.Subchunk2ID[1] = 'd';
  // outputChunk.Subchunk2ID[2] = 'a';
  // outputChunk.Subchunk2ID[3] = 't';
  // outputChunk.Subchunk2ID[4] = 'a';
