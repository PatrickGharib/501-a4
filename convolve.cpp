#include <ctime>
#include <array>
#include <iostream>
#include <string>
#include <fstream>
#include <vector> 
#include <cstdint>
#include <algorithm>
#include <string.h>
using namespace std;



//struct to organize wave info from input 
typedef struct WAVEHEADER{
    //TODO change to uint32_t and stuff for these cause its apparently really bad 
    
    uint8_t RIFF[4];
    uint32_t ChunkSize;
    uint8_t format[4];
    uint8_t Subchunk1ID[4];
    uint32_t Subchunk1Size;
    uint16_t AudioFormat;
    uint16_t NumChannels;
    uint32_t SampleRate; 
    uint32_t ByteRate;
    uint16_t BlockAlign;
    uint16_t BitsPerSample;
    uint8_t Subchunk2ID[4];
    uint32_t Subchunk2Size;
} wavHdr;

//Prototyping
void convolve(double x[], int N, double h[], int M, double y[], int P);

int getWaveFileSize(FILE *inFile);
int8_t* getWavData(FILE *wavFile);

int main(int argc, char** argv){

    //check user input for right input
    if(argc != 4){
        cout << "Usage: convolve <inputfile(.wav)> <IRfile(.wav)> <outputFile(.wav)>"; 
        exit(0);
    }

    wavHdr wavHeader1;
    wavHdr wavHeader2;

    int sizeOfHeader1 = sizeof(wavHeader1), lengthOfFile = 0;
    int sizeOfHeader2 = sizeof(wavHeader2);
    
    const char* inputFile = argv[1];


    //TODO need to move this to save some mem so that is not just doubleing
    FILE* wavFile = fopen(inputFile, "r");

    if (wavFile == NULL){
        fprintf(stderr, "Unable to open inputFile: %s\n", inputFile);
        return 1;
    }

    //read in the header
    size_t bytesRead = fread(&wavHeader1, 1, sizeOfHeader1, wavFile);
    uint16_t bytesPerSample1 = wavHeader1.BitsPerSample/8;
    int64_t numberOfSamples1 = wavHeader1.ChunkSize/bytesPerSample1;
    //if (bytesRead > 0){
    static const uint64_t BUFFERSIZE = wavHeader1.Subchunk2Size;
    int8_t* buffer1 = new int8_t[BUFFERSIZE];
    while((bytesRead = fread(buffer1, sizeof buffer1[0], BUFFERSIZE / (sizeof buffer1[0]), wavFile)) > 0){}
    double* wavData1 = new double[BUFFERSIZE];
        //TODO make this into a method as a refactoring
    for(int i = 0; i < BUFFERSIZE; i++){wavData1[i] = buffer1[i];}
        //delete [] buffer1;
        //buffer1 = nullptr;
    //fileLength = getFileSize(wavFile);
    //} 



   const char* iRFile = argv[2];

   wavFile = fopen(iRFile, "r");

    if (wavFile == NULL){
        fprintf(stderr, "Unable to open IRFile: %s\n", iRFile);
        return 1;
    }
    //read in the header
    bytesRead = fread(&wavHeader2, 1, sizeOfHeader2, wavFile);
    uint16_t bytesPerSample2 = wavHeader2.BitsPerSample/8;
    int64_t numberOfSamples2 = wavHeader2.ChunkSize/bytesPerSample2;
   // if (bytesRead > 0){
        //TODO extract this into a method so that 
        static const uint64_t BUFFERSIZE2 = wavHeader2.Subchunk2Size;
        int8_t* buffer2 = new int8_t[BUFFERSIZE2];
        //TODO add vector record the data part of wave
        while((bytesRead = fread(buffer2, sizeof buffer2[0], BUFFERSIZE2 / (sizeof buffer2[0]), wavFile)) > 0){}
       // delete [] buffer1;
       // buffer1 = nullptr;
        //fileLength = getFileSize(wavFile);
      //}
    fclose(wavFile);
    
     
    int64_t outputNumberOfSamples = numberOfSamples1 + numberOfSamples2 - 1;
    double* outputArray = new double[outputNumberOfSamples];
    
    double* wavData2 = new double[BUFFERSIZE2];
    for(int i = 0; i< BUFFERSIZE2;i++){wavData2[i] = buffer2[i];}
    convolve(wavData1, numberOfSamples1, wavData2, numberOfSamples2, outputArray, outputNumberOfSamples);
    
    wavHdr outputWaveHeader;
    outputWaveHeader.RIFF[0] =  'R';
    outputWaveHeader.RIFF[1] =  'I';
    outputWaveHeader.RIFF[2] =  'F';
    outputWaveHeader.RIFF[3] =  'F';

    outputWaveHeader.format[0] = 'W';
    outputWaveHeader.format[1] = 'A';
    outputWaveHeader.format[2] = 'V';
    outputWaveHeader.format[3] = 'E';

    outputWaveHeader.Subchunk1ID[0] = 'f';
    outputWaveHeader.Subchunk1ID[1] = 'm';
    outputWaveHeader.Subchunk1ID[2] = 't';

    outputWaveHeader.Subchunk2ID[1] = 'd';
    outputWaveHeader.Subchunk2ID[2] = 'a';
    outputWaveHeader.Subchunk2ID[3] = 't';
    outputWaveHeader.Subchunk2ID[4] = 'a';

    outputWaveHeader.Subchunk2Size = wavHeader1.NumChannels * outputNumberOfSamples * (wavHeader1.BitsPerSample/8);
    outputWaveHeader.Subchunk1Size = 16;
    outputWaveHeader.AudioFormat = 1;
    outputWaveHeader.NumChannels = wavHeader1.NumChannels;
    outputWaveHeader.ChunkSize = 36 + outputWaveHeader.Subchunk2Size;
    outputWaveHeader.BlockAlign = wavHeader1.NumChannels * (wavHeader1.BitsPerSample/8);
    outputWaveHeader.ByteRate = wavHeader1.SampleRate * outputWaveHeader.BlockAlign;
    outputWaveHeader.BitsPerSample = wavHeader1.BitsPerSample;
    outputWaveHeader.SampleRate = wavHeader1.SampleRate;
    
    cout << outputArray;
    int sizeOfHeader3 = sizeof(outputWaveHeader);
    const char* outputFile = argv[3];
    wavFile = fopen(outputFile, "w"); 
    fwrite(&outputWaveHeader,1, sizeOfHeader3,wavFile);
    fwrite(outputArray,outputWaveHeader.BitsPerSample/8,outputWaveHeader.Subchunk2Size/(outputWaveHeader.BitsPerSample/8),wavFile);
    fclose(wavFile);
   
    
    
    
    return 0;

}

int getWaveFileSize(FILE *inFile){
    int fileSize = 0; 
    fseek(inFile, 0, SEEK_END);
    fileSize = ftell(inFile); fseek(inFile,0, SEEK_SET); 
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

void convolve(double x[], int N, double h[], int M, double y[], int P)
{
  int n, m;
   
  /*  Make sure the output buffer is the right size: P = N + M - 1  */
  if (P != (N + M - 1)) {
    printf("Output signal vector is the wrong size\n");
    printf("It is %-d, but should be %-d\n", P, (N + M - 1));
    printf("Aborting convolution\n");
    
    return;
  }
  /*  Clear the output buffer y[] to all zero values  */  
  for (n = 0; n < P; n++){
    y[n] = 0.0;
    //cout << y[n];
    }
  int start;
  int stop;
  bool flag = true;  
  /*  Do the convolution  */
  /*  Outer loop:  process each input value x[n] in turn  */
  for (n = 0; n < N; n++) {
    if (flag){
      start = clock();
      flag = false;
    }
    /*  Inner loop:  process x[n] with each sample of h[]  */
    for (m = 0; m < M; m++)
      {y[n+m] += x[n] * h[m];
      
      //cout << y[n+m];
     }

     //calculate estimated time to completion
     if (n % 10000 == 0){
      stop = clock();
      cout<<"estimated time : " << (((N-n)/10000)*((stop - start)/double(CLOCKS_PER_SEC)))/60 << endl;
      flag = true;
     }
  }
}
