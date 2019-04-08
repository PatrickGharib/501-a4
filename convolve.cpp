

#include <iostream>
#include <string>
#include <fstream>
#include <vector> 

using namespace std;

//struct to organize wave info from input 
typedef struct WAVEHEADER{
    //TODO change to uint32_t and stuff for these cause its apparently really bad 
    
    uint8_t RIFF[4];
    uint32_t ChunkSize;
    char format[4];
    char Subchunk1ID;
    uint32_t Subchunk1Size;
    uint16_t AudioFormat;
    uint16_t NumChannels;
    uint16_t SampleRate; 
    uint16_t ByteRate;
    uint16_t BlockAlign;
    uint16_t BitsPerSample;
    uint8_t Subchunk2ID[4];
    uint16_t Subchunk2Size;
}wavHdr;
int getWaveFileSize(FILE* fileSize);
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
    //const char* iRFile = argv[2];im
    FILE* wavFile = fopen(inputFile, "r");
    if (wavFile = nullptr){
        fprintf(stderr, "Unable to open inputFile: %s\n", filePath);
        return 1;
    }

    size_t bytesRead fread(&fread(&wavHeader,1, sizeOfheader, waveFile));
    if (bytesRead > 0){
        static const uint64_t BUFFERSIZE = 4096;
        int8_t buffer = new int8_t[BUFFERSIZE];
        vector wavData
        while((bytesRead = fread(buffer, sizeof buffer[0], BUFFERSIZE/ (sizeof buffer[0]), wavFile)) > 0){
            
        }
    }

    string inputFile = argv[1];
    string IRfile = argv[2];
    string outputFile = argv[3];
    cout << inputFile << " " << IRfile << " " << outputFile;



}
int getWaveFileSize(File* fileSize);
void getWaveData(){

}