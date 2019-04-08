#include <iostream>
#include <string>
#include <fstream>
#include <vector> 

using namespace std;

//struct to organize wave info from input 
typedef struct WAVEHEADER{
    //TODO change to uint32_t and stuff for these cause its apparently really bad 
    
    char RIFF[4];
    unsigned int ChunkSize;
    char format[4];
    char Subchunk1ID; 
    unsigned long Subchunk1Size;
    unsigned short AudioFormat;
    unsigned short NumChannels;
    unsigned long SampleRate; 
    unsigned long ByteRate;
    unsigned short BlockAlign;
    unsigned long BitsPerSample;
    char Subchunk2ID[4];
    unsigned long Subchunk2Size;
}wavHdr;
int main(int argc, char** argv){
    //check user input for right input
    if(argc != 4){
        cout << "Usage: convolve <inputfile(.wav)> <IRfile(.wav)> <outputFile(.wav)>"; 
        exit(0);
    }

    wavhdr wavHead;
    int sizeOfheader = sizeof(wavHeader), lengthOfFile = 0;
    const char* inputFile = argv[1];
    //TODO need to move this to save some mem so that is not just floating
    //const char* iRFile = argv[2];
    FILE* wavFile = fopen

    string inputFile = argv[1];
    string IRfile = argv[2];
    string outputFile = argv[3];
    cout << inputFile << " " << IRfile << " " << outputFile;



}
int getWaveFileSize(File* fileSize);
void getWaveData(){

}